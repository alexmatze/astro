#!/usr/bin/python
#====================================================
#This script is used to successfully extract gaussian 
#components in a radio map. 
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

path_darks_b="darks_before.txt"
list_darks_b=list(open(path_darks_b,"r"))
n_darks_b=len(list_darks_b)

path_flats_b="flats_before.txt"
list_flats_b=list(open(path_flats_b,"r"))
n_flats_b=len(list_flats_b)

path_darks_a="darks_after.txt"
list_darks_a=list(open(path_darks_a,"r"))
n_darks_a=len(list_darks_a)

path_flats_a="flats_after.txt"
list_flats_a=list(open(path_flats_a,"r"))
n_flats_a=len(list_flats_a)


path_data="data.txt"
list_data=list(open(path_data,"r"))
n_data=len(list_data)

# Hier noch die Positionen und Radien eintragen
target_pos=[1,1] #[y,x]
target_rad=[20]
cali1_pos=[1,1]
cali1_rad=[12]
cali2_pos=[1,1]
cali2_rad=[12]
cali3_pos=[1,1]
cali3_rad=[12]

dither=50 #um wie viel das maximum von bild zu bild variieren kann

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

darks_data_b=[]
darks_head_b=[]
darks_data_a=[]
darks_head_a=[]

flats_data_b=[]
flats_head_b=[]
flats_data_a=[]
flats_head_a=[]


for i in range(n_darks_b):
	path=list_darks_b[i]
	path=path.rstrip()
	hdulist=pyfits.open(path) * 1
	dataset=hdulist[0].data 
	header=hdulist[0].header 
	darks_data_b.append(dataset)
	darks_head_b.append(header)





