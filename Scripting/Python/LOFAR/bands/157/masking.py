#!/usr/bin/python
#====================================================
#This script is used to identify regions
#in fits images by using modified SEED+GROW algo
#SEED and GROW procedure iterate until satisfied
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

path_in="band4.fits"
path_out="mask.fits"
#path_resid="band7_resid2.fits"
path_resid="modeltest.fits"
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

def seed(image,mask,increment,stopval): #Find Seeds to grow where Seed is found as isolated maximum within max and increment but >stopval
	tmp_mask=mask * 1	#create a temporary mask for region/noregion assoc
	canlist=[]		#open empty list for seed candidates
	for y in range(mask.shape[0]):	#fill temporary mask with 0:any region and 1:no region
		for x in range(mask.shape[1]):
			if mask[y,x]==0:
				tmp_mask[y,x]=1
			else:
				tmp_mask[y,x]=0
	tmp_image= image * tmp_mask #create temporary image with blanked out region pixels
	maxval=np.max(tmp_image)	#identify remaining maximum
	if maxval <= stopval:	#if no maximum above limit is found -> no new seed
		print "No further Seeds to plant"
		return mask,maxval
        for y in range(mask.shape[0]):	#identify candidate pixels with max-inc<val<max
                for x in range(mask.shape[1]):
			if tmp_image[y,x]>maxval-maxval*increment:
				canlist.append([y,x])
	tmplst=[]
	for candidate in canlist:
		if candidate[0]>1 and candidate[1]>1 and candidate[0]<mask.shape[0]-1 and candidate[1]<mask.shape[1]-1: #drop edge pixel candidates
			if mask[candidate[0]+1,candidate[1]]==0 and mask[candidate[0],candidate[1]+1]==0 and mask[candidate[0]-1,candidate[1]]==0 and mask[candidate[0],candidate[1]-1]==0:
				tmplst.append(candidate) #only keep candidates not adjacent to regions
	canlist=tmplst * 1
	reg_no=np.max(mask)+1
	seedlist=[]
	#collect connected candidates to regions
	flr_mask = mask * 0.
	for can in canlist:
		flr_mask[can[0],can[1]]=-1
	for y in range(mask.shape[0]):	#search for candidade pixel and fill adjecient candidates with next no
		for x in range(mask.shape[1]):
			if flr_mask[y,x]==-1:
				flr_mask[y,x]=reg_no * 1
				if flr_mask[y+1,x] == -1 or flr_mask[y-1,x]==-1 or flr_mask[y,x+1]==-1 or flr_mask[y,x-1]==-1:
					condition = True
					chck=[[y,x]]
					new_chck=[]
					tsl=[]
					i=0
					while condition:
						for cent in chck:
							if flr_mask[cent[0]+1,cent[1]]==-1:
								new_chck.append([cent[0]+1,cent[1]])
								flr_mask[cent[0]+1,cent[1]]=reg_no * 1
								tsl.append([cent[0]+1,cent[1]+0])
							if flr_mask[cent[0]-1,cent[1]]==-1:
								new_chck.append([cent[0]-1,cent[1]])
								flr_mask[cent[0]-1,cent[1]]=reg_no * 1
								tsl.append([cent[0]-1,cent[1]+0])
							if flr_mask[cent[0],cent[1]+1]==-1:
								new_chck.append([cent[0],cent[1]+1])
								flr_mask[cent[0],cent[1]+1]=reg_no * 1
								tsl.append([cent[0]+0,cent[1]+1])
							if flr_mask[cent[0],cent[1]-1]==-1:
								new_chck.append([cent[0],cent[1]-1])
								flr_mask[cent[0],cent[1]-1]=reg_no * 1
								tsl.append([cent[0]+0,cent[1]-1])

						chck = new_chck *1
						new_chck = []
						if len(new_chck) == 0:
							seedlist.append(tsl * 1)
							condition = False

				else:
					flr_mask[y,x]=reg_no * 1
					seedlist.append([[y,x]])
				
				reg_no = reg_no +1

	#only consider maximum within each region
	for seeds in seedlist:
		maxseed=0.
		maxx=0
		maxy=0
		for seed in seeds:
			if image[seed[0],seed[1]] > maxseed:
				maxy=seed[0]
				maxx=seed[1]
				maxseed = image[seed[0],seed[1]]
		mask[maxy,maxx]=np.max(mask)+1
	
	#write new mask with additional seeds and return mask
	return mask,maxval

def grow(image,mask,maxval_p,increment,stopval):
	#do growing procedure
	grow_condition=True
	new_mask=mask*1
	loopcnt=0
	maxval=maxval_p * 1.
	while grow_condition:
		old_mask = mask * 1
		for y in range(mask.shape[0]):	#scan for existing region pixels and let them grow within parameters
			for x in range(mask.shape[1]):	
				if mask[y,x]!=0:
					if y>0 and x>0 and y<mask.shape[0] and x<mask.shape[1]:
						if mask[y-1,x]==0 and image[y-1,x]>stopval and image[y-1,x]>maxval*(1-increment):
							mask[y-1,x]=mask[y,x]
						if mask[y+1,x]==0 and image[y+1,x]>stopval and image[y+1,x]>maxval*(1-increment):
							mask[y+1,x]=mask[y,x]
						if mask[y,x-1]==0 and image[y,x-1]>stopval and image[y,x-1]>maxval*(1-increment):
							mask[y,x-1]=mask[y,x]
						if mask[y,x+1]==0 and image[y,x+1]>stopval and image[y,x+1]>maxval*(1-increment):
							mask[y,x+1]=mask[y,x]
		print "Finished Scan!"
		if old_mask.all()==mask.all():
			blend_mask=mask*1
			for y in range(mask.shape[0]):
				for x in range(mask.shape[1]):
					if mask[y,x]!=0.:
						blend_mask[y,x]=0.
					else:
						blend_mask[y,x]=1.
			maxvaln=np.max(blend_mask*image)
			grow_condition=False


	return mask,maxvaln
			
	

#=====================================================
#                      Main
#=====================================================



#read in all your file:

hdulist=pyfits.open(path_in) * 1
dataset=hdulist[0].data[0][0] * 1
header=hdulist[0].header


hdulist_resid=pyfits.open(path_resid) * 1
dataset_resid=hdulist_resid[0].data[0][0] * 1

low_level = 0.005
grow_loops = 20
init_mask=dataset * 0.
g_mask = dataset * 0. + 1.
loop_count = 0
loop_cond = True
run_mask=init_mask * 1


while loop_cond:
	print loop_count
	print "Seeds get planted"
	new_mask,maximumval = seed(dataset,run_mask,0.99,low_level)
	print "Seeding complete!"
	print "Start growing process in "+str(grow_loops)+" iterations!"
	for i in range(grow_loops):
		print "Grwoth round "+str(i)
		new_mask,maximumval = grow(dataset,new_mask,maximumval,0.3,low_level)
	if new_mask.all() == run_mask.all():
		final_mask = new_mask
		loop_cond = False
	loop_count = loop_count + 1
		


#two-region fix
reg_no_split = 93

tmp_mask = final_mask * 1
topfl=0.
topflr=0.
botfl=0.
botflr=0.

for y in range(tmp_mask.shape[0]):
	for x in range(tmp_mask.shape[1]):
		if tmp_mask[y,x]==0.:
			final_mask[y,x]=0

		elif tmp_mask[y,x]>0. and tmp_mask[y,x]<=reg_no_split:
			final_mask[y,x]=1
			botfl=botfl+dataset[y,x]
			botflr=botflr+dataset_resid[y,x]

		elif tmp_mask[y,x]>0. and tmp_mask[y,x]>reg_no_split:
			final_mask[y,x]=2
			topfl=topfl+dataset[y,x]
			topflr=topflr+dataset_resid[y,x]




print "Measured Flux values:"
print "--------------- FOV:"
print "F_fov: "+str(np.sum(dataset))+"Jy"
print "--------------- Original:"
print "Top: "+str(topfl)+"Jy"
print "Bot: "+str(botfl)+"Jy"

print "--------------- Residual:"
print "Top: "+str(topflr)+"Jy"
print "Bot: "+str(botflr)+"Jy"

print "--------------- Ratio HS+CHS/Core:"
print str((topflr+botflr)/(topfl-topflr))

print "--------------- Ratio Core/HS+CHS:"
print str((topfl-topflr)/(topflr+botflr))





if os.path.isfile(path_out):
	os.system("rm " + path_out)
	
pyfits.writeto(path_out, final_mask , header)










