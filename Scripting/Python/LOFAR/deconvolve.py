import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from dateutil import parser
from scipy import fftpack
import os
import astropy.io.fits as pyfits

#######################################
#############USER-INTERFACE############
#######################################
input_path_image='0836+710_resid2.fits'
output_path_image='0836+710_resid2_decon.fits'




#######################################
############PROGRAM-CORE###############
#######################################

def convolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft*psf_fft)))

def deconvolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/psf_fft)))

def ddeconvolve(star, psf, epsilon):
    starfft= fftn(star)
    psffft= fftn(psf) + epsilon
    deconvolved = ifftn(starfft/psffft)
    deconvolved = np.abs(deconvolved)
    return(deconvolved)

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def set0tomin64(img):
    for y in range(len(img)):
        for x in range(len(img[0])):
            if img[y,x] < np.finfo(np.float32).tiny:
                img[y,x] = np.finfo(np.float32).tiny
    return img
	
def setinfto064(img):
    for y in range(len(img)):
        for x in range(len(img[0])):
            if np.isinf(img[y,x]):
                img[y,x] = 0.
    return img
	
def setnanto064(img):
    for y in range(len(img)):
        for x in range(len(img[0])):
            if np.isnan(img[y,x]):
                img[y,x] = 0.
    return img


#Data-readin
#readin image data
hdul = fits.open(input_path_image)
neu = hdul[0].data
image = neu #[0][0] #Jy

#readin of the beam
bmaj = hdul[0].header['BMAJ'] #deg
bmin = hdul[0].header['BMIN'] #deg
bpa = hdul[0].header['BPA'] #deg
bpar = bpa * np.pi /180.

#readin of y-x-Sizes
x_size=hdul[0].header["NAXIS1"]
y_size=hdul[0].header["NAXIS2"]

center_x=int(np.ceil(x_size/2.0))
center_y=int(np.ceil(y_size/2.0))

#readin of the axis scale
scale = hdul[0].header['CDELT1'] #deg/pixel

#readin of the difmap noise
auto_noise = hdul[0].header['NOISE']

#readin of the observation date
date = parser.parse(hdul[0].header['DATE-OBS'])
date = date.strftime('%Y-%m-%d')

#readin default Source name as saved in fitsfile and set sourcename
name_default = hdul[0].header['OBJECT']

#build image of beam with same imagesize as image
beam_img=image * 0.
for x in range(x_size):
	for y in range(y_size):
		beam_img[y,x]=twoD_Gaussian((x,y),1,center_x,center_y,bmin/scale, bmaj/scale, -bpar,0)
beam_img = beam_img / beam_img.sum()

#fix arrays to np-def
image=set0tomin64(np.array(image,dtype=np.float64))
beam_img=set0tomin64(np.array(beam_img,dtype=np.float64))



#plt.figure("Beam"); plt.imshow(beam_img,norm=colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
#plt.show()



if os.path.isfile(output_path_image):
	os.system("rm " + output_path_image)

pyfits.writeto(output_path_image,beam_img)

image_decon=ddeconvolve(image,beam_img,0)

real_image_decon = np.real(image_decon)
real_image_decon=setinfto064(real_image_decon)
real_image_decon=setnanto064(real_image_decon)
real_image_decon=set0tomin64(real_image_decon)

plt.figure("Beam"); plt.imshow(real_image_decon,norm=colors.LogNorm());plt.colorbar();plt.gca().invert_yaxis()
plt.show()
