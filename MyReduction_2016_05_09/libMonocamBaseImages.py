# -*- coding: utf-8 -*-
"""
Created on Thu May 26 10:56:23 2016

@author: dagoret-campagnesylvie
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy import units as u
import ccdproc

import os

NB_OF_CHANNELS=16


#-------------------------------------------------------------------------------------
def BuildFilelist(path,name,ext='.fits',start=1,stop=99):
    '''
    Make the list of filenames required by ccdproc
    
    input:
       path : path of files
       name : common root of bias filenames
       ext  : extension of filenames
       start,stop : indexes of files
    output:
       full filename list
    '''
    filelist = []
    for num in range(start,stop+1,1):
        strnum=biasnumberstr= '{0:02d}'.format(num)  # python >= 2.6
        filename=name+strnum+ext
        fullfilename=os.path.join(path,filename)
        filelist.append(fullfilename)
    return filelist
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
def SaveCCDListIntoFitsFile(ccdlist,outfitsfilename,meta,imagetyp='master_bias'):
	'''
	SaveCCDListIntoFitsFile::
		
		inputs : 
		- ccdlist : list of CCDData images to be written in the file
		- outfitsfilename : string of the output fits filename
		- meta : the information we want to keep in the header
		output:
		  outfitsfilename : output fits file			
	'''
	# first create the HDU list and the primary  HDU (header)
	
	
	hdul = fits.HDUList() 
	hdul.append(fits.PrimaryHDU()) 

	header=hdul[0].header

	header=meta
	header['IMAGETYP']=imagetyp
	#header['DATE-ANA']=date : this info should already be in the input header

	print 'header to be written in file ::'
	print '-------------------------------'
	header

	index=0 # channel index
	# loop on CCD channels
	for ccdchan in ccdlist :
	    index=index+1
	    hd=ccdchan.header
	    print index
	    print hd
	    img=ccdchan.data
	    hdul.append(fits.ImageHDU(data=img)) 

	hdul.info()
	hdul.writeto(outfitsfilename,clobber=True) 

     #hdul.close()

#----------------------------------------------------------------------------------------------------
def oscan_and_trim(image_list):
    """
    Remove overscan and trim a list of images. The original list is replaced by a list of images
    with the changes applied.

    Implementation done by ccdproc


    Parameters:
    ----------

    image_list :: List of CCDData corresponding to images
    
    
    """
    for idx, img in enumerate(image_list):
        oscan = ccdproc.subtract_overscan(img,overscan=img[:,521:544], add_keyword={'oscan_sub': True, 'calstat': 'O'}, model=models.Polynomial1D(1))
        image_list[idx] = ccdproc.trim_image(oscan[:,10:521], add_keyword={'trimmed': True, 'calstat': 'OT'})
#------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------
def bn_median(masked_array, axis=None):
    """
    Perform fast median on masked array
    
    Parameters
    ----------
    
    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.
    
    axis : int, optional
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
    """
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))
					
#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
def avg_over_images(masked_arr, axis=0):
    """
    Calculate average pixel value along specified axis
    """
    return ma.mean(masked_arr, axis=axis)
#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
def med_over_images(masked_arr, axis=0):
    """
    Calculate median pixel value along specified axis
    
    Uses bottleneck.nanmedian for speed
    """
    
    dat = masked_arr.data.copy()
    dat[masked_arr.mask] = np.NaN
    return bn.nanmedian(dat, axis=axis)
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
def overscan_trim_and_sigma_clip_median(image_list, clip_baseline_func=med_over_images):
    """
    Combine a list of images using median
    
    This function does several steps:
    
    1. Subtract overscan
    2. Trim image
    3. sigma clip image using a median of the unclipped stack as the baseline
    4. combine the images on the list using median
    
    ** It modifies the images in the input list. **
    """
    oscan_and_trim(image_list)
    combo = ccdproc.Combiner(image_list)
    combo.sigma_clipping(func=clip_baseline_func)
    return combo
#-------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------
def ShowImagesSet(ccdlist):
    '''
    Shows the whole set of CCD images
     - inputs argument:
       path : path of the fits file
       filename of the fits file
     - output the images of the whole CCD   
    '''
     
    NX=8 # number of images along the horizontal axis
    NY=2 # number of images along the vertical axis
    
    f, axarr = plt.subplots(NY,NX,sharex='col', sharey='row',figsize=(15,15)) # figure organisation
    
    f.subplots_adjust(hspace=0.125,wspace=0.1)

    for index in range(NB_OF_CHANNELS):  
        ix=index%8
        iy=index/8
        image_data = ccdlist[index].data
        V_MIN=image_data.flatten().min()
        V_MAX=image_data.flatten().max()
        im=axarr[iy,ix].imshow(image_data,vmin=V_MIN,vmax=V_MAX)  # plot the image
        if ix==0 and iy==0:
            im0=im
        plottitle='channel {}'.format(index+1)
        axarr[iy,ix].set_title(plottitle)
    
    title='Master Biases'
    cax = f.add_axes([0.95, 0.12, 0.03, 0.78]) # [left,bottom,width,height]    
    f.colorbar(im0, cax=cax)
   
    plt.suptitle(title,size=16)
    plt.savefig('viewccd.pdf', bbox_inches='tight')
#----------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------
#   Main to test the library					
#-------------------------------------------------------------------------------------------------
if __name__ == "__main__":
	in_masterbias_filename='masterbias1.fits'
	out_masterbias_filename='masterbias1_outtest.fits'
	
	hdu_list = fits.open(in_masterbias_filename)
	hdu_list.info()

	header=hdu_list[0].header
	meta=header
	
	print '- Header read from file :: '
	print header
	
	# all CCDPROC data collector : each channel as a list of biases data
	allccd = []
	for chan in range(1,NB_OF_CHANNELS+1,1):
		ccd_chan =  ccdproc.CCDData(hdu_list[chan].data,meta=header,unit="adu")
		allccd.append(ccd_chan)

	SaveCCDListIntoFitsFile(allccd,out_masterbias_filename,meta)
	hdu_list.close()
