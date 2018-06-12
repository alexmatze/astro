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

path_darks_b="darks_b.txt"
list_darks_b=list(open(path_darks_b,"r"))
n_darks_b=len(list_darks_b)

path_flats_b="flats_b.txt"
list_flats_b=list(open(path_flats_b,"r"))
n_flats_b=len(list_flats_b)

path_darks_a="darks_a.txt"
list_darks_a=list(open(path_darks_a,"r"))
n_darks_a=len(list_darks_a)

path_flats_a="flats_a.txt"
list_flats_a=list(open(path_flats_a,"r"))
n_flats_a=len(list_flats_a)


path_data="data.txt"
list_data=list(open(path_data,"r"))
n_data=len(list_data)


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
                print("ERROR: No Data in Array!!!")
                return False
        y_len=np.shape(a[0])[0]
        x_len=np.shape(a[0])[1]

        for i in range(len(a)):
                test= (np.shape(a[0]) == np.shape(a[i]))
                if test == False:
                        print("ERROR: Sizes are different!!!")
                        return False

        c=np.zeros((y_len,x_len))
        for x in range(x_len):
                for y in range(y_len):
                        vals=[]
                        for i in range(len(a)):
                                vals.append(a[i][y,x])
                        c[y,x]= np.median(vals) * 1.

        return c

def triangular (r1, r2, r3, x, y,*print_check):
        for i in r1:
                i = i * 1.
        for i in r2:
                i = i * 1.
        for i in r3:
                i = i * 1.
        x = x * 1.
        y = y * 1.
#       if print_check: print(r1,r2,r3,x,y)
        a = ( ((x-r1[2])*(r3[1]-r1[1])-(y-r1[1])*(r3[2]-r1[2]))/((r2[2]-r1[2])*(r3[1]-r1[1])-(r2[1]-r1[1])*(r3[2]-r1[2])) ) *1.
        b = ( ((x-r1[2])*(r2[1]-r1[1])-(y-r1[1])*(r2[2]-r1[2]))/((r3[2]-r1[2])*(r2[1]-r1[1])-(r3[1]-r1[1])*(r2[2]-r1[2])) ) *1.


        val = ( r1[0] + a*(r2[0]-r1[0]) + b*(r3[0]-r1[0]) ) *1.


        return val




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
        dataset=hdulist[0].data * 1
        header=hdulist[0].header
        darks_data_b.append(dataset)
        darks_head_b.append(header)

for i in range(n_darks_a):
        path=list_darks_a[i]
        path=path.rstrip()
        hdulist=pyfits.open(path) * 1
        dataset=hdulist[0].data * 1
        header=hdulist[0].header
        darks_data_a.append(dataset)
        darks_head_a.append(header)

for i in range(n_flats_b):
        path=list_flats_b[i]
        path=path.rstrip()
        hdulist=pyfits.open(path) * 1
        dataset=hdulist[0].data * 1
        header=hdulist[0].header
        flats_data_b.append(dataset)
        flats_head_b.append(header)

for i in range(n_flats_a):
        path=list_flats_a[i]
        path=path.rstrip()
        hdulist=pyfits.open(path) * 1
        dataset=hdulist[0].data * 1
        header=hdulist[0].header
        flats_data_a.append(dataset)
        flats_head_a.append(header)

# We start to calculate median flats for before/after measurement

x_size=darks_head_b[0]["NAXIS1"]
y_size=darks_head_b[0]["NAXIS2"]
exp_time_darks=darks_head_b[0]["EXPTIME"]
exp_time_flats=flats_head_b[0]["EXPTIME"]

if os.path.isfile("./geterrorimage/median_darks_b.fits"):
        print("Median flats and darks already exist and are read in.")

        hdulist=pyfits.open("./geterrorimage/median_darks_b.fits") * 1
        median_darks_b=hdulist[0].data * 1

        hdulist=pyfits.open("./geterrorimage/median_darks_a.fits") * 1
        median_darks_a=hdulist[0].data * 1

        hdulist=pyfits.open("./geterrorimage/median_flats_b.fits") * 1
        median_flats_b=hdulist[0].data * 1

        hdulist=pyfits.open("./geterrorimage/median_flats_a.fits") * 1
        median_flats_a=hdulist[0].data * 1

else:
        print("Medians flats and darks do not yet exist. Calculating Median flats and darks.")
        median_darks_b=median(darks_data_b)/exp_time_darks *1
        median_darks_a=median(darks_data_a)/exp_time_darks *1

        median_flats_b=median(flats_data_b)/exp_time_flats *1
        median_flats_a=median(flats_data_a)/exp_time_flats *1

        median_flats_b=(median_flats_b - median_darks_b) *1
        median_flats_a=(median_flats_a - median_darks_a) *1

        new_head_d=darks_head_b[0]
        new_head_f=flats_head_b[0]
        new_head_d["EXPTIME"]= 1.0
        new_head_f["EXPTIME"]= 1.0
        pyfits.writeto("./geterrorimage/median_darks_b.fits",median_darks_b,new_head_d)
        pyfits.writeto("./geterrorimage/median_darks_a.fits",median_darks_a,new_head_d)
        pyfits.writeto("./geterrorimage/median_flats_b.fits",median_flats_b,new_head_f)
        pyfits.writeto("./geterrorimage/median_flats_a.fits",median_flats_a,new_head_f)
        print("Medians flats and darks were calculated and saved.")

# Compute Error Images for Median Flats and Darks

median_darks_b_error = np.sqrt(median_darks_b) *1
median_darks_a_error = np.sqrt(median_darks_a) *1
median_flats_b_error = np.sqrt(median_flats_b) *1
median_flats_a_error = np.sqrt(median_flats_a) *1


# now let's correct each target-image 

time=[]
date=[]
mjd=[]


# we need to know, when last image was taken:

path=list_data[n_data-1]
path=path.rstrip()
hdulist=pyfits.open(path) * 1
header=hdulist[0].header
final_mjd=Time([header["DATE-OBS"]+"T"+header["TIME-OBS"]],format="isot",scale="utc").mjd

box = 201  #Boxgroesse ungerade waehlen!
radius = int((box - 1) * 0.5)
over = 20

for i in range(n_data):
        print(i)

        path=list_data[i]
        path=path.rstrip()
        print("Current image:")
        print(path)
        corr_path="./geterrorimage/cor_"+path
        error_path="./geterrorimage/error_"+path
        bias_path="./geterrorimage/bias_"+path
        hdulist=pyfits.open(path) * 1
        dataset=hdulist[0].data *1
        #header=hdulist[0].header *1
        header=hdulist[0].header
        exp_time=header["EXPTIME"]

        # Define Error of the original dataset
        dataset_error = np.sqrt(dataset) * 1
        dataset=dataset/exp_time * 1
        dataset_error = dataset_error/exp_time * 1

        time.append(header["TIME-OBS"])
        date.append(header["DATE-OBS"])
        time_mjd=Time([header["DATE-OBS"]+"T"+header["TIME-OBS"]],format="isot",scale="utc").mjd
        mjd.append(time_mjd)

        if i==0:
                first_mjd=time_mjd

        # Define corrected darks and flats as well as their errors
        corr_darks= median_darks_b + (median_darks_a - median_darks_b) * (time_mjd-first_mjd)/(final_mjd-first_mjd)
        corr_darks_error = np.sqrt((median_darks_b_error * (1 - (time_mjd-first_mjd)/(final_mjd-first_mjd)))**2 + (median_darks_a_error * ((time_mjd-first_mjd)/(final_mjd-first_mjd)))**2)
        corr_flats= median_flats_b + (median_flats_a - median_flats_b) * (time_mjd-first_mjd)/(final_mjd-first_mjd)
        corr_flats_error = np.sqrt((median_flats_b_error * (1 - (time_mjd-first_mjd)/(final_mjd-first_mjd)))**2 + (median_flats_a_error * ((time_mjd-first_mjd)/(final_mjd-first_mjd)))**2)

        # Correct the dataset for darks and flats and propagate the error
        dataset = (dataset - corr_darks) / corr_flats
        dataset_error = (1/corr_flats)*np.sqrt(corr_darks_error**2 + dataset_error**2 + (dataset*corr_flats_error)**2)
        dataset = extrapolate(dataset, 1654, 160)
        dataset = extrapolate(dataset, 681, 3476)


        if os.path.isfile(bias_path):
                if i==0:
                        print("Bias images already exist and are read in.")

                hdulist=pyfits.open(bias_path) * 1
                bias_image=hdulist[0].data * 1

        else:
                if i==0:
                        print("Bias images do not exist yet. Calculating bias images (this can take a while..)")


                #------------------------------
                # Correct Image for the bias
                #------------------------------
                bias_image = dataset * 1
                bias_image = bias_image * 0.


                # Boxgroesse und Overlap definieren

                #box = 201  #Boxgroesse ungerade waehlen!
                #radius = (box - 1) * 0.5
                #over = 20
                checkdist = box * 3.
                #x_size = [4096]
                #y_size =[4096]
                x_size=header["NAXIS1"]
                y_size=header["NAXIS2"]
                pos = [radius , radius]
                median_value=[]


                # Anzahl der Boxen berechnen:

                nx_box = int(((x_size - box)/(box - over) + 1)) *1
                ny_box = int(((y_size - box)/(box - over) + 1)) *1



                # For loop ueber alle Boxen - Mediane der Boxen berechnen und den Mittelpunkt auf diesen Wert setzen:
                for j in range(ny_box):
                        for i in range(nx_box):

                                # Define district, compute median and set center pixel to that median value:
                                district = dataset[pos[0]-radius:pos[0]+radius,pos[1]-radius:pos[1]+radius] * 1
                                temp_median = np.median(district)
                                median_value.append([temp_median * 1., pos[0] * 1, pos[1] *1])





                                # Shift center pixel position
                                pos[1] = ( pos[1] + box - over ) * 1

                        # Define additional box for the right side of the image:
                        pos[1] = ( x_size - radius ) * 1
                        district = dataset[pos[0]-radius:pos[0]+radius,pos[1]-radius:pos[1]+radius] * 1
                        temp_median = np.median(district)
                        median_value.append([temp_median * 1., pos[0] * 1, pos[1] *1])





                        # Shift center pixel position
                        pos[1] = radius
                        pos[0] = ( pos[0] + box - over ) * 1

                #Define additional boxes for the top side of the image:


                pos[0] = (y_size - radius ) * 1

                # For loop for the bottom row:
                for k in range(nx_box):
                        district = dataset[pos[0]-radius:pos[0]+radius,pos[1]-radius:pos[1]+radius] * 1
                        temp_median = np.median(district)
                        median_value.append([temp_median * 1., pos[0] * 1, pos[1] *1])





                        # Shift center pixel position
                        pos[1] = ( pos[1] + box - over ) * 1


                # Define additional box for the right side of the image:


                pos[1] = ( x_size - radius ) * 1
                district = dataset[pos[0]-radius:pos[0]+radius,pos[1]-radius:pos[1]+radius] * 1
                temp_median = np.median(district)
                median_value.append([temp_median * 1., pos[0] * 1, pos[1] *1])

                # Set the value of the middle pixel:

                for val, y, x in median_value:

                        bias_image[y,x] = val


                # Find points with largest x/y value

                for indx, (val, y, x) in enumerate(median_value):
                        if indx==0:
                                max_x = int(x) * 1
                        elif x > max_x:
                                max_x = int(x) * 1

                for indx, (val, y, x) in enumerate(median_value):
                        if indx==0:
                                max_y = int(y) * 1
                        elif x > max_y:
                                max_y = int(y) * 1


                # Find points with smallest x/y value

                for indx, (val, y, x) in enumerate(median_value):
                        if indx==0:
                                min_x = int(x) *1
                        elif x < min_x:
                                min_x = int(x) *1

                for indx, (val, y, x) in enumerate(median_value):
                        if indx==0:
                                min_y = int(y) *1
                        elif x < min_y:
                                min_y = int(y) *1

                # Vorbereitung fuer das Fuellen des bias-image

                condition=True
                current_point = [0,min_y,min_x]
                secondpoint=[0,0,0]
                thirdpoint=[0,0,0]
                fourthpoint=[0,0,0]
                secondpoint[1]=current_point[1]
                secondpoint[2]=(current_point[2]+box-over)
                thirdpoint[1]=(current_point[1]+box-over)
                thirdpoint[2]=current_point[2]
                fourthpoint[1]=(current_point[1]+box-over)
                fourthpoint[2]=(current_point[2]+box-over)
                secondp=[0,0,0]
                thirdp=[0,0,0]
                fourthp=[0,0,0]

                #While-Schleife fuellt das bias-image:

                while condition:
                        secondpoint[1]=current_point[1]
                        secondpoint[2]=(current_point[2]+box-over)
                        thirdpoint[1]=(current_point[1]+box-over)
                        thirdpoint[2]=current_point[2]
                        fourthpoint[1]=(current_point[1]+box-over)
                        fourthpoint[2]=(current_point[2]+box-over)

                        if (current_point[2]+box-over)>max_x and (current_point[1]+box-over)< max_y:

                                secondp[2]= max_x *1
                                secondp[1]= current_point[1] *1
                                thirdp[2]= current_point[2] *1
                                thirdp[1]= (current_point[1]+box-over) *1
                                fourthp[2]= max_x *1
                                fourthp[1]= (current_point[1]+box-over) *1

                                for indx2, (val2, y2, x2) in enumerate(median_value):
                                        if current_point[1]==y2 and current_point[2]==x2:
                                                current_point[0]=val2 *1.
                                        elif secondp[1]==y2 and secondp[2]==x2:
                                                secondp[0]=val2 *1.
                                        elif thirdp[1]==y2 and thirdp[2]==x2:
                                                thirdp[0]=val2 *1.
                                        elif fourthp[1]==y2 and fourthp[2]==x2:
                                                fourthp[0]=val2 *1.
                                y_list = range(current_point[1]+1, thirdp[1]+1)
                                x_list = range(current_point[2]+1, secondp[2]+1)
                                m_stg = (1.*thirdp[1]-1.*current_point[1])/(1.*secondp[2]-1.*current_point[2])
                                for y in y_list:
                                        for x in x_list:
                                                m_stg_p = (1.*y-1.*current_point[1])/(1.*x-1.*current_point[2])
                                                if m_stg_p > m_stg :
                                                        bias_image[y,x] = triangular(current_point,thirdp,fourthp,x,y)
                                                else:
                                                        bias_image[y,x] = triangular(current_point,secondp,fourthp,x,y)

                                current_point[1]= (current_point[1]+box-over) *1
                                current_point[2]= min_x


                        elif (current_point[1]+box-over)> max_y and (current_point[2]+box-over)<max_x:

                                secondp[2]= (current_point[2]+box-over) *1
                                secondp[1]= current_point[1] *1
                                thirdp[2]= current_point[2] *1
                                thirdp[1]= max_y *1
                                fourthp[2]= (current_point[2]+box-over) *1
                                fourthp[1]= max_y *1

                                for indx2, (val2, y2, x2) in enumerate(median_value):
                                        if current_point[1]==y2 and current_point[2]==x2:
                                                current_point[0]=val2 *1.
                                        elif secondp[1]==y2 and secondp[2]==x2:
                                                secondp[0]=val2 *1.
                                        elif thirdp[1]==y2 and thirdp[2]==x2:
                                                thirdp[0]=val2 *1.
                                        elif fourthp[1]==y2 and fourthp[2]==x2:
                                                fourthp[0]=val2 *1.
                                y_list = range(current_point[1]+1, thirdp[1]+1)
                                x_list = range(current_point[2]+1, secondp[2]+1)
                                for y in y_list:
                                        for x in x_list:
                                                if ((y-current_point[1])/(thirdp[1]-current_point[1])) > ((x-current_point[2])/(secondp[2]-current_point[2])):
                                                        bias_image[y,x] = triangular(current_point,thirdp,fourthp,x,y)
                                                else:
                                                        bias_image[y,x] = triangular(current_point,secondp,fourthp,x,y)
                                current_point[2] = (current_point[2] + box - over) * 1


                        elif(current_point[1]+box-over)> max_y and (current_point[2]+box-over)>max_x:

                                secondp[2]= max_x *1
                                secondp[1]= current_point[1] *1
                                thirdp[2]= current_point[2] *1
                                thirdp[1]= max_y *1
                                fourthp[2]= max_x *1
                                fourthp[1]= max_y *1

                                for indx2, (val2, y2, x2) in enumerate(median_value):
                                        if current_point[1]==y2 and current_point[2]==x2:
                                                current_point[0]=val2 *1.
                                        elif secondp[1]==y2 and secondp[2]==x2:
                                                secondp[0]=val2 *1.
                                        elif thirdp[1]==y2 and thirdp[2]==x2:
                                                thirdp[0]=val2 *1.
                                        elif fourthp[1]==y2 and fourthp[2]==x2:
                                                fourthp[0]=val2 *1.
                                y_list = range(current_point[1]+1, thirdp[1]+1)
                                x_list = range(current_point[2]+1, secondp[2]+1)
                                for y in y_list:
                                        for x in x_list:
                                                if ((y-current_point[1])/(thirdp[1]-current_point[1])) > ((x-current_point[2])/(secondp[2]-current_point[2])):
                                                        bias_image[y,x] = triangular(current_point,thirdp,fourthp,x,y)
                                                else:
                                                        bias_image[y,x] = triangular(current_point,secondp,fourthp,x,y)
                                condition=False


                        else:
                                for indx2, (val2, y2, x2) in enumerate(median_value):
                                        if current_point[1]==y2 and current_point[2]==x2:
                                                current_point[0]=val2 *1.
                                        elif secondpoint[1]==y2 and secondpoint[2]==x2:
                                                secondpoint[0]=val2 *1.
                                        elif thirdpoint[1]==y2 and thirdpoint[2]==x2:
                                                thirdpoint[0]=val2 *1.
                                        elif fourthpoint[1]==y2 and fourthpoint[2]==x2:
                                                fourthpoint[0]=val2 *1.

                                y_list = range(current_point[1]+1, thirdpoint[1]+1)
                                x_list = range(current_point[2]+1, secondpoint[2]+1)
                                for y in y_list:
                                        for x in x_list:
                                                if (y-current_point[1]) > (x-current_point[2]):
                                                        bias_image[y,x] = triangular(current_point,thirdpoint,fourthpoint,x,y,False)
                                                else:
                                                        bias_image[y,x] = triangular(current_point,secondpoint,fourthpoint,x,y,False)
                                current_point[2] = (current_point[2] + box-over) * 1

                if os.path.isfile(bias_path):
                        os.system("rm " + bias_path)
                new_head_b=darks_head_b[0]
                pyfits.writeto(bias_path,bias_image,new_head_b)

        # Substract bias from the data
        dataset = dataset - bias_image

        # Aeussere Pixel abschneiden, da das bias_image nur das Zentrum modelliert
        x3list = range(0,x_size)
        y3list = range(0,y_size)
        for x3 in x3list:
                for y3 in y3list:
                        if x3 < radius + 1 or y3 < radius + 1  or x3 > x_size - radius - 1 or y3 > y_size - radius - 1:
                                dataset[y3,x3] = dataset[y3,x3]*0.
                                dataset_error[y3,x3] = dataset_error[y3,x3]*0.




        if os.path.isfile(corr_path):
                os.system("rm " + corr_path)
        pyfits.writeto(corr_path, dataset, header)

        if os.path.isfile(error_path):
                os.system("rm " + error_path)
        pyfits.writeto(error_path, dataset_error, header)
