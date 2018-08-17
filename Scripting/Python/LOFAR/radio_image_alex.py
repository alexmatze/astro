import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from dateutil import parser
#######################################
#############USER-INTERFACE############
#######################################
band_name="4"
freq_name="126-160"
#freq_name="126-160_Gaussfit"
#freq_name="126-160_stack"



#input path of the interferometric FITS-file, and the model-fitsfile and output path of the resulting image
#input_path_image = '/scratch/local/akappes/git/astro/Scripting/Python/LOFAR/bands/'+freq_name+'/band'+band_name+'_resid2.fits'
input_path_image = '/scratch/local/akappes/git/astro/Scripting/Python/LOFAR/0836+710_resid2.fits'
#input_path_image = '/scratch/local/akappes/git/astro/Scripting/Python/LOFAR/0836+710_comp2.fits'
#input_path_image = '/scratch/local/akappes/git/astro/Scripting/Python/LOFAR/0836+710_stacked_hba.fits'
input_path_model = input_path_image
output_path = './'+freq_name+'.pdf'

#Source name
#name = '0836+710 @'+freq_name+' MHz (Band '+band_name+')'
name = '' #'0836+710 @'+freq_name+' MHz'
name_color = 'Grey' # color in which the source-name and the date is written

#noise and lowest level-cut (sigma)
noise = 0 #[JY] if noise <= 0, the noise for the hole map from difmap is used
sigma = 3.1 # if sigma <= 0, 3 is used

#unit-selection,map limits and axe ratio
unit = 'arcsec' #possible units: 'mas', 'arcsec', 'arcmin' and 'deg'. default (any other string) value is 'deg'. Limits are controlled in this unit.
ra_min = -2
ra_max = 1.5
dec_min = -2.5
dec_max =1.5

#image look tuning
#contour plot
contour = True # if True a contour plot is done
contour_color = ['Grey'] #input: array of color-strings if None the contour-colormap (contour_cmap) is used
contour_cmap = None #matplotlib colormap string
contour_alpha = 1 # Transparency
contour_width = 0.5 # contour linewidth
#Range of Color/grey scale
scale_min = -0.0398212879 #if 'None': defined automatically
scale_max = 1.4932113886 #if 'None': defined automatically

#image colormap
im_colormap =True # if True a image colormap is done
im_color = 'viridis' #string for matplotlib colormap

#model overplot
#overplot Gaussian
overplot_gauss = False #if True all Gaussian components are plottet
gauss_linewidth = 0.5 # linewidth of the ellipse
gauss_color = ['Black'] # Array of the colors of the single ellipses, if the number of color entries is smaler than the numbers of ellipses, every ellipse has the color of the first entry

#overplot clean
overplot_clean = False #if True all clean components are plottet
clean_size = None # float to set the sympol size of the clean components; None sets default value
clean_color = 'Grey' # string for sympol color
clean_alpha = 1 # float for sympol transparency
clean_linewidth = 0.5 #clean linewidth of the symbol


#######################################
############PROGRAM-CORE###############
#######################################


#Data-readin
#readin image data
hdul = fits.open(input_path_image)
neu = hdul[0].data
if len(neu)==1:
	image = neu[0][0]
else:
	image = neu



#readin of the beam
bmaj = hdul[0].header['BMAJ'] #deg
bmin = hdul[0].header['BMIN'] #deg
bpa = hdul[0].header['BPA'] #deg

#readin of the axis scale
scale = hdul[0].header['CDELT1'] #deg/pixel

#readin of the difmap noise
auto_noise = hdul[0].header['NOISE']

#readin of the observation date
date = parser.parse(hdul[0].header['DATE-OBS'])
date = ""#date.strftime('%Y-%m-%d')

#readin default Source name as saved in fitsfile and set sourcename
name_default = hdul[0].header['OBJECT']
if name == None:
	name = name_default

# #readin of model components
# hdul2 = fits.open(input_path_model)
# mod = hdul2[1].data
# mod = np.array(mod)
# mod = mod.view((mod.dtype[0], len(mod.dtype.names)))
# mod = np.transpose(mod, [1,0])

# #sorting of Gaussian and Clean
# cond_g = mod[3] >0
# cond_c = mod[3] == 0

# g_x = np.extract(cond_g, mod[1])
# g_y = np.extract(cond_g, mod[2])
# g_maj = np.extract(cond_g, mod[3])
# g_min = np.extract(cond_g, mod[4])
# g_pos = np.extract(cond_g, mod[5])

# c_x = np.extract(cond_c, mod[1])
# c_y = np.extract(cond_c, mod[2])


#Adjustments
#level adjustment
if noise <= 0:
	noise = auto_noise

if sigma <= 0:
	sigma = 3

level0 = noise*sigma

#unit selection and adjustment
if unit == 'arcmin':
	scale = scale*60
	bmaj = bmaj*60
	bmin = bmin*60
	# g_x = g_x*60
	# g_y = g_y*60
	# g_maj = g_maj*60
	# g_min = g_min*60
	# c_x = c_x*60
	# c_y = c_y*60

if unit == 'arcsec':
	scale = scale*60*60
	bmaj = bmaj*60*60
	bmin = bmin*60*60
	# g_x = g_x*60*60
	# g_y = g_y*60*60
	# g_maj = g_maj*60*60
	# g_min = g_min*60*60
	# c_x = c_x*60*60
	# c_y = c_y*60*60

if unit == 'mas':
	scale = scale*60*60*1000
	bmaj = bmaj*60*60*1000
	bmin = bmin*60*60*1000
	# g_x = g_x*60*60*1000
	# g_y = g_y*60*60*1000
	# g_maj = g_maj*60*60*1000
	# g_min = g_min*60*60*1000
	# c_x = c_x*60*60*1000
	# c_y = c_y*60*60*1000

x = np.linspace(-len(image[0])*0.5*scale,(len(image[0])*0.5-1)*scale,len(image[0]))
y = np.linspace(len(image[0])*0.5*scale,-(len(image[0])*0.5-1)*scale,len(image[0]))

extent = np.max(x), np.min(x), np.min(y), np.max(y)
fig, ax = plt.subplots()

#PLOTTING
#plot look tuning
axe_ratio = 'scaled'
plt.axis(axe_ratio)
plt.xlim(ra_min,ra_max)
plt.ylim(dec_min,dec_max)
plt.gca().invert_xaxis()
plt.xlabel('Relative RA ['+unit+']')
plt.ylabel('Relative DEC ['+unit+']')

if scale_min == 'None':
	scale_min = -level0
if scale_max == 'None':
	scale_max = 0.5*np.max(image)
print("Applied scale range from "+str(scale_min)+" to "+str(scale_max))


#image colormap
if im_colormap == True:
	col = ax.imshow(image, cmap=im_color,norm=colors.SymLogNorm(linthresh=level0, linscale=0.5,vmin=scale_min, vmax=scale_max),extent=extent,origin='lower')
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = fig.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy]')

#contour plot
if contour == True:
	ax.contour(x,y,image,linewidths=contour_width,levels=[-level0,level0,level0*2,level0*2**2,level0*2**3,level0*2**4,level0*2**5,level0*2**6,level0*2**7,level0*2**8,level0*2**9,level0*2**10], colors=contour_color, alpha=contour_alpha,cmap=contour_cmap, norm=colors.SymLogNorm(linthresh=level0, linscale=0.5,vmin=scale_min, vmax=scale_max))
print("Defined contour-levels:"+str([-level0,level0,level0*2,level0*2**2,level0*2**3,level0*2**4,level0*2**5,level0*2**6,level0*2**7,level0*2**8,level0*2**9,level0*2**10]))
# #overplot
# #overplot clean
# if overplot_clean == True:
# 	ax.scatter(c_x,c_y,marker='_', c=clean_color, s=clean_size, alpha=clean_alpha, lw =clean_linewidth, zorder=2)
# 	ax.scatter(c_x,c_y,marker='|', c=clean_color, s=clean_size, alpha=clean_alpha, lw =clean_linewidth, zorder=2)

# #overplot gaus
# if len(gauss_color) < len(g_x):
# 	gauss_color = np.full(len(g_x),gauss_color[0])

# if overplot_gauss == True:
# 	for i in range(len(g_x)):
# 		#plotting ellipses
# 		e = Ellipse([g_x[i],g_y[i]],g_maj[i],g_min[i],-g_pos[i]+90, fill=False, zorder=2,color=gauss_color[i],lw=gauss_linewidth)
# 		ax.add_artist(e)

# 		#plotting crosses of the ellipses
# 		maj1_x = g_x[i]-np.sin(-np.pi/180*g_pos[i])*g_maj[i]*0.5
# 		maj1_y = g_y[i]+np.cos(-np.pi/180*g_pos[i])*g_maj[i]*0.5
# 		maj2_x = g_x[i]+np.sin(-np.pi/180*g_pos[i])*g_maj[i]*0.5
# 		maj2_y = g_y[i]-np.cos(-np.pi/180*g_pos[i])*g_maj[i]*0.5

# 		min1_x = g_x[i]-np.sin(-np.pi/180*(g_pos[i]+90))*g_min[i]*0.5
# 		min1_y = g_y[i]+np.cos(-np.pi/180*(g_pos[i]+90))*g_min[i]*0.5
# 		min2_x = g_x[i]+np.sin(-np.pi/180*(g_pos[i]+90))*g_min[i]*0.5
# 		min2_y = g_y[i]-np.cos(-np.pi/180*(g_pos[i]+90))*g_min[i]*0.5


# 		ax.plot([maj1_x,maj2_x],[maj1_y,maj2_y], color = gauss_color[i], lw = gauss_linewidth)
# 		ax.plot([min1_x,min2_x],[min1_y,min2_y], color = gauss_color[i], lw = gauss_linewidth)


#set achse label default
if (unit != 'mas') and (unit != 'arcsec') and (unit != 'arcmin'):
	unit = 'deg'




#set beam ellipse, Sourcename and observation date possitions
size_x = np.absolute(ra_max)+np.absolute(ra_min)
size_y = np.absolute(dec_max)+np.absolute(dec_min)
ell_dist = 1
if size_x>size_y:
	ell_x = ra_max-bmaj*ell_dist
	ell_y = dec_min+bmaj*ell_dist
	name_x = ra_max-size_x*0.05
	name_y = dec_max-size_x*0.05
	date_x = ra_min+size_x*0.05
	date_y = name_y
else:
	ell_x = ra_max-bmaj*ell_dist
	ell_y = dec_min+bmaj*ell_dist
	name_x = ra_max-size_y*0.05
	name_y = dec_max-size_y*0.05
	date_x = ra_min+size_y*0.05
	date_y = name_y

#plot beam-ellipse
e = Ellipse([ell_x,ell_y],bmaj,bmin,-bpa+90, fc='grey',zorder=2)
ax.add_artist(e)

#plot date
ax.text(date_x,date_y,date,color=name_color,ha='right',va='top')

#plot name
ax.text(name_x,name_y,name,color=name_color,ha='left',va='top')

# tight layout
plt.tight_layout()

#print(g_pos)
# save image
plt.savefig(output_path,bbox_inches='tight')
