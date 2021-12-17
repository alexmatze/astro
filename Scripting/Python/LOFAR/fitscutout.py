#Fitsfile Imageclipper based on the work by E. Bonnassieux
#https://github.com/ebonnassieux/Scripts/blob/master/FitsCutoutFromFits.py
#Python3 adoption, dealing with different fits image formats, rotated fields

import sys
import argparse
import numpy as np
from astropy.io import fits
from astLib.astWCS import WCS
from astLib.astCoords import calcRADecSearchBox
from astLib.astImages import clipRotatedImageSectionWCS
from astLib.astImages import saveFITS

def readArguments():
    # create parser
    parser=argparse.ArgumentParser(description="Make a FITS cutout from a FITS image.")
    parser.add_argument("--filename",type=str,help="Name of the .fits file you want a cutout of",required=True)
    parser.add_argument("--outname",type=str,help="Name of the output .fits file, optional",required=False)#,nargs=argparse.REMAINDER)
    parser.add_argument("--RA",metavar="HH:MM:SS",type=str, help="Right Ascension in hexadecimal hour angle",required=True)
    parser.add_argument("--Dec",metavar="HH:MM:SS",type=str,help="Declination in hexadecimal degrees",required=True)
    parser.add_argument("--size",metavar="arcmin",type=float,default=5,help="Size of cutout, in arcmin. Default is 5")
    # parse
    args=parser.parse_args()
    return vars(args)

def MakeCutout(filename,RA,dec,ArcMinSize,cutoutname=None):
    # parse RA, Dec coordinates
    if type(RA)==str:
        print("Cropped image centre:")
        print("RA   : %s"%RA)
        print("Dec  : %s"%dec)
        print("Size : %s '"%ArcMinSize)
        RAdeg  = HHMMSStoDegrees(RA)*15. # this converts RA from hours to degrees
        Decdeg = HHMMSStoDegrees(dec)
    imhdu = fits.open(filename)
    imwcs = WCS(filename)
    if (imhdu[0].data).ndim ==4:
        imdata = imhdu[0].data[0,0,:,:]
    elif (imhdu[0].data).ndim ==2:
        imdata = imhdu[0].data
    else:
        print("Your input FITS-file does not match the required shape. The datacube should have a dimension of 2 or 4.")
        return 1
    # make cutout box
    rmin,rmax,dmin,dmax=calcRADecSearchBox(RAdeg,Decdeg,ArcMinSize/60.)
    cutout = clipRotatedImageSectionWCS(imdata,imwcs,RAdeg,Decdeg,2.*ArcMinSize/60.)
    im=cutout["data"]
    if cutoutname==None:
        cutoutname=filename+".cutout.fits"
    saveFITS(cutoutname,cutout['data'],cutout['wcs'])
    print("Cutout is: %s"%cutoutname)
    imhdu.close()


def HHMMSStoDegrees(HHMMSS):
   # convert HHMMSS string to angular value float
   HH,MM,SS=np.array(HHMMSS.split(":")).astype(float)
   degrees=HH+MM/60.+SS/3600.
   return degrees

if __name__=="__main__":
    # parse arguments for this function
    args=readArguments()
    # assign variables
    filename=args["filename"]
    outname=args["outname"]
    ra=args["RA"]
    dec=args["Dec"]
    arcminsize=args["size"]
    # launch script
    MakeCutout(filename,ra,dec,arcminsize,outname)
