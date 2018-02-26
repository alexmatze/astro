#!/usr/local/bin/python

import pyfits as cunt
import numpy as np
import os

CRVAL3
names = []
baender=[]

namen = ['band2','band3','band4','band5','band7','band8','band9','band10','band11','band12','band13']

calib_val=[0.403701,0.367172,0.414706,0.307139,0.260366,0.157054,0.140854,0.174828,0.350102,0.348418,0.377940]

for i in range(2,6):
	names.append('/scratch/local/LC4/26/BAND'+str(i)+'_IMAGE/BAND'+str(i)+'_REMAP.fits')

for i in range(7,14):
	names.append('/scratch/local/LC4/26/BAND'+str(i)+'_IMAGE/BAND'+str(i)+'_REMAP.fits')


for path in names:
	baender.append(cunt.open(path))
	
for i in range(len(baender)):
	baender[i][0].data = baender[i][0].data * calib_val[i]
	baender[i][1].data['FLUX'] = baender[i][1].data['FLUX'] * calib_val[i]
	outname = './corrected/'+namen[i]+'.fits'
	if os.path.isfile(outname):
		os.system('rm '+outname)
	cunt.writeto('./corrected/'+namen[i]+'.fits',baender[i][0].data,baender[i][0].header)
	cunt.append('./corrected/'+namen[i]+'.fits',baender[i][1].data,baender[i][1].header)



outm = (baender[0][0].data[0][0] * 1) * 0.0
outm_median = (baender[0][0].data[0][0] * 1) * 0.0

for i in range(256):
	for j in range(256):
		dummy=[]
		for obj in baender:
			outm[i][j] = outm[i][j] + obj[0].data[0][0][i][j]
			dummy.append(obj[0].data[0][0][i][j])
		outm_median[i][j] = np.median(dummy)

outm = outm / (1.0 * len(baender))
outm_med_mean=outm-outm_median

outname1="testout.fits"
outname2="testout_median.fits"
outname3="testout_res.fits"

if os.path.isfile(outname1):
	os.system("rm "+outname1)
if os.path.isfile(outname2):
	os.system("rm "+outname2)
if os.path.isfile(outname3):
	os.system("rm "+outname3)

cunt.writeto('testout.fits',outm,baender[3][0].header)
cunt.writeto('testout_median.fits',outm_median,baender[3][0].header)
cunt.writeto('testout_res.fits',outm_med_mean,baender[3][0].header)

#----------------------------------- Adding here the section for Spectral map 


