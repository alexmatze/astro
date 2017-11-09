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
	
for i in range(n_darks_a):
	path=list_darks_a[i]
	path=path.rstrip()
	hdulist=pyfits.open(path) * 1
	dataset=hdulist[0].data 
	header=hdulist[0].header 
	darks_data_a.append(dataset)
	darks_head_a.append(header)
	
for i in range(n_flats_b):
	path=list_flats_b[i]
	path=path.rstrip()
	hdulist=pyfits.open(path) * 1
	dataset=hdulist[0].data 
	header=hdulist[0].header 
	flats_data_b.append(dataset)
	flats_head_b.append(header)
		
for i in range(n_flats_a):
	path=list_flats_a[i]
	path=path.rstrip()
	hdulist=pyfits.open(path) * 1
	dataset=hdulist[0].data 
	header=hdulist[0].header 
	flats_data_a.append(dataset)
	flats_head_a.append(header)
	
# We start to calculate median flats for before/after measurement

x_size=darks_head_b[0]["NAXIS1"]
y_size=darks_head_b[0]["NAXIS2"]
exp_time_darks=darks_head_b[0]["EXPTIME"]
exp_time_flats=flats_head_b[0]["EXPTIME"]

median_darks_b=median(darks_data_b)
median_darks_a=median(darks_data_a)

median_darks_b=median_darks_b
median_darks_a=median_darks_a

median_flats_b=median(flats_data_b)
median_flats_a=median(flats_data_a)

median_flats_b=median_flats_b/exp_time_flats - median_darks_b/exp_time_darks
median_flats_a=median_flats_a/exp_time_flats - median_darks_a/exp_time_darks

# now let's correct each target-image and calculate the Flux

time=[]
date=[]
mjd=[]
target_flux=[]
cali1_flux=[]
cali2_flux=[]
cali3_flux=[]
calibrated_flux=[]

# we need to know, when last image was taken:

path=list_data[n_data-1]
path=path.rstrip()
hdulist=pyfits.open(path) * 1
header=hdulist[0].header
final_mjd=Time([header["DATE-OBS"]+"T"+header["TIME-OBS"]],format="isot",scale="utc").mjd

for i in range(n_data):

		
	target_pos=[1,1] #[y,x] Correct Positions for the Image
	cali1_pos=[1,1]
	cali2_pos=[1,1]
	cali3_pos=[1,1]
	
	path=list_data[i]
	path=path.rstrip()
	corr_path="./corrected/cor_"+path
	hdulist=pyfits.open(path) * 1
	dataset=hdulist[0].data 
	header=hdulist[0].header
	exp_time=header["EXPTIME"]
	
	time.append(header["TIME-OBS"])
	date.append(header["DATE-OBS"])
	time_mjd=Time([header["DATE-OBS"]+"T"+header["TIME-OBS"]],format="isot",scale="utc").mjd
	mjd.append(time_mjd)

	if i==0:
		first_mjd=time_mjd

	corr_darks= median_darks_b + (median_darks_a - median_darks_b) * (time_mjd-first_mjd)/(final_mjd-first_mjd)
	corr_flats= median_flats_b + (median_flats_a - median_flats_b) * (time_mjd-first_mjd)/(final_mjd-first_mjd)

	dataset = (dataset/exp_time - corr_darks/exp_time_darks) #/ corr_flats
	#dataset = extrapolate(dataset, 1654, 160)
	#dataset = extrapolate(dataset, 681, 3476)
	
	#correct each image for the background noise
	dataset = dataset - np.median(dataset)
	
	if os.path.isfile(corr_path):
		os.system("rm " + corr_path)
	
	pyfits.writeto(corr_path, dataset, header)
	if i==0:
		pos_max=np.where(dataset == np.max(dataset))
	else:
		area=dataset[pos_max[0]-dither:pos_max[0]+dither,pos_max[1]-dither:pos_max[1]+dither] * 1
		new_pos=np.where(area==np.max(area))
		new_pos_a = pos_max[0] - (dither + 1 -new_pos[0]) 
		new_pos_b = pos_max[1] - (dither + 1 -new_pos[1]) 
		pos_max=(new_pos_a,new_pos_b)*1

	print "dataset: ",
	print path
	print "x: ",
	print pos_max[1],
	print "		y: ",
	print pos_max[0]

	target_pos=[pos_max[0] * 1,pos_max[1] * 1] #[y,x] Correct Positions for the Image
	cali1_pos=[target_pos[0] + 164 ,target_pos[1] - 90]
	cali2_pos=[target_pos[0] + 130 ,target_pos[1] + 460]
	cali3_pos=[target_pos[0] - 226 ,target_pos[1] + 128]
	

	#here calc for flux stuff!
	target_index=[]
	for j in range((target_pos[1]-target_rad[0]),(target_pos[1]+target_rad[0]+1)):
		for k in range((target_pos[0]-target_rad[0]),(target_pos[0]+target_rad[0]+1)):
			if dist(target_pos[1],target_pos[0],j,k) <= target_rad[0]:
				target_index.append([k,j])	

	cali1_index=[]
	for j in range((cali1_pos[1]-cali1_rad[0]),(cali1_pos[1]+cali1_rad[0]+1)):
		for k in range((cali1_pos[0]-cali1_rad[0]),(cali1_pos[0]+cali1_rad[0]+1)):
			if dist(cali1_pos[1],cali1_pos[0],j,k) <= cali1_rad[0]:
				cali1_index.append([k,j])	
	cali2_index=[]
	for j in range((cali2_pos[1]-cali2_rad[0]),(cali2_pos[1]+cali2_rad[0]+1)):
		for k in range((cali2_pos[0]-cali2_rad[0]),(cali2_pos[0]+cali2_rad[0]+1)):
			if dist(cali2_pos[1],cali2_pos[0],j,k) <= cali2_rad[0]:
				cali2_index.append([k,j])	

	cali3_index=[]
	for j in range((cali3_pos[1]-cali3_rad[0]),(cali3_pos[1]+cali3_rad[0]+1)):
		for k in range((cali3_pos[0]-cali3_rad[0]),(cali3_pos[0]+cali3_rad[0]+1)):
			if dist(cali3_pos[1],cali3_pos[0],j,k) <= cali3_rad[0]:
				cali3_index.append([k,j])	

	tmp_flux=0.
	for j in target_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]	
	target_flux.append(tmp_flux * 1.)

	tmp_flux=0.
	for j in cali1_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]	
	cali1_flux.append(tmp_flux * 1.)


	tmp_flux=0.
	for j in cali2_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]	
	cali2_flux.append(tmp_flux * 1.)


	tmp_flux=0.
	for j in cali2_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]	
	cali3_flux.append(tmp_flux * 1.)
		
	
for i in range(len(target_flux)):
	calibrated_flux.append(target_flux[i]/(cali1_flux[i]+cali2_flux[i]+cali3_flux[i]))	

data_output=np.array(zip(time,date,target_flux,cali1_flux,cali2_flux,cali3_flux,calibrated_flux), dtype=[("time","S16"), ("date","S16"), ("target_flux",float), ("cali1_flux",float), ("cali2_flux",float), ("cali3_flux",float),("calibrated_flux",float)])

outname="fluxes_data.txt"
np.savetxt(outname, data_output, fmt=["%s"]*2 + ["%.5f"]*5, delimiter="|")
print "Your data is saved in '%s'" % (outname)









