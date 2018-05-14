import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from dateutil import parser
from scipy import fftpack

#######################################
#############USER-INTERFACE############
#######################################
input_path_image=''
output_path_image=''




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

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel() / g.sum()


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
beam_img=twoD_Gaussian(image,1,center_x,center_y,bmin/scale, bmaj/scale, bpar,0)


image_decon=deconvolve(image,beam_img)
