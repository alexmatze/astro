#!/usr/bin/env python
import os,sys
import glob
import pyrap.tables as pt
import numpy as np

########################################################################
def main(ms_input, SkymodelPath, Radius="5.", DoDownload="True"):
    """
    Download the TGSS skymodel for the target field
    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    SkymodelPath : str
        Full name (with path) to the skymodel; if YES is true, the TGSS skymodel will be downloaded here
    Radius : string with float (default = "5.")
        Radius for the TGSS cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Download or not the TGSS skymodel.
        "Force": download skymodel from TGSS, delete existing skymodel if needed.
        "True" or "Yes": use existing skymodel file if it exists, download skymodel from
                         TGSS if it does not.
        "False" or "No": Do not download skymodel, raise an exception if skymodel
                         file does not exist.
    """
  FileExists = os.path.isfile(SkymodelPath)
    if (not FileExists and os.path.exists(SkymodelPath)):
        raise ValueError("download_tgss_skymodel_target: WTF! Path: \"%s\" exists but is not a file!"%(SkymodelPath))
    download_flag = False
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(SkymodelPath)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print "USING the exising skymodel in "+ SkymodelPath
            return
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print "USING the exising skymodel in "+ SkymodelPath
            return
         else:
            raise ValueError("download_tgss_skymodel_target: Path: \"%s\" does not exist and TGSS download is disabled!"%(SkymodelPath))

    # If we got here, then we are supposed to download the skymodel.
    assert download_flag == True # Jaja, belts and suspenders...
    print "DOWNLOADING TGSS Skymodel for the target into "+ SkymodelPath

    # Reading a MS to find the coordinate (pyrap)
    [RATar,DECTar]=grab_coord_MS(input2strlist_nomapfile(ms_input)[0])

    # Downloading the skymodel
    os.system("wget -O "+SkymodelPath+ " \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv2.cgi?coord="+str(RATar)+","+str(DECTar)+"&radius="+Radius+"&unit=deg&deconv=y\' ")

    if not os.path.isfile(SkymodelPath):
        raise IOError("download_tgss_skymodel_target: Path: \"%s\" does not exist after trying to download TGSS skymodel."%(SkymodelPath))

    return


########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=' Download the TGSS skymodel for the target field')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which a TGSS skymodel will be download.')
    parser.add_argument('SkyTar', type=str,
                        help='Full name (with path) to the skymodel; the TGSS skymodel will be downloaded here')
    parser.add_argument('--Radius', type=float,
                        help='Radius for the TGSS cone search in degrees')

    args = parser.parse_args()
    radius=5
    if args.Radius:
        radius=args.Radius

    main(args.MSfile,args.SkyTar, str(radius))
