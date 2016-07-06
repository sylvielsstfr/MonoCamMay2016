
# coding: utf-8

# # Extraction of One spectrum
# =============================================
# 
# - Creation : Tuesday 2016 June 14th
# - Update / July 5th 2016
# - Author Sylvie Dagoret-Campagne 
# - affiliation : LAL/IN2P3/CNRS
# 
# Study spectrum of HD158485
# 
########################################################
# ## 1) Import package



import matplotlib.pyplot as plt
import numpy as np

from astropy.modeling import models
from astropy import units as u
from astropy import nddata
from astropy.io import fits

import ccdproc
print 'ccdproc version',ccdproc.__version__

from astropy.modeling import models

import photutils
from astropy.stats import sigma_clipped_stats
from photutils import daofind
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D

from scipy import stats 
from scipy import ndimage
import os
from datetime import datetime, timedelta

from scipy.interpolate import UnivariateSpline

import sys


import libMonocamBaseImages           # my tool library written to do that CCD reduction

#------------------------------------------------------------------------------
def SeparateSpectra(inspectra,x0):
    '''
    Cut the two spectra
    '''
    rightspectra=inspectra[x0:]
    revleftspectra=inspectra[:x0]
    leftspectra=   revleftspectra[::-1]
    #rightspectra=rightspectra[np.where(rightspectra>0)]
    #leftspectra=leftspectra[np.where(leftspectra>0)]
    
    return leftspectra,rightspectra
#-----------------------------------------------------------------------------    

#----------------------------------------------------------------------------- 
def GetSpectrumBackground(inspectra,start,stop,skybg):
    '''
    Return the background    
    '''
    cropedbg=inspectra[start:stop]
    purebg=cropedbg[np.where(cropedbg!=skybg)]  # remove region of the bing star
    
    return purebg
#----------------------------------------------------------------------------- 

#----------------------------------------------------------------------------- 
def DiffSpectra(spec1,spec2,bg):
    '''
    Make the difference of the tow spectra 
    
    '''
    N1=spec1.shape[0]
    N2=spec2.shape[0]
    N=np.min([N1,N2])
    spec1_croped=spec1[0:N]
    spec2_croped=spec2[0:N]
    diff_spec=np.average((spec1_croped-spec2_croped)**2)/bg**2
    return diff_spec  
#----------------------------------------------------------------------------- 
 
#---------------------------------------------------------------------------   
def FindCenter(fullspectrum,xmin,xmax,specbg):
    '''
    This allows to find the minimum of the chi2 to get the center
    '''
    
    all_x0=np.arange(xmin,xmax,1)
    NBPOINTS=np.shape(all_x0)
    chi2=np.zeros(NBPOINTS)
    for idx,x0 in np.ndenumerate(all_x0):
        spec1,spec2=SeparateSpectra(fullspectrum,x0)
        chi2[idx]=DiffSpectra(spec1,spec2,specbg)
    return all_x0,chi2  
#---------------------------------------------------------------------------      

#--------------------------------------------------------------------------- 
def DiffAmplitudes(sp1,sp2,gain,basebg):
    '''
    get the differences of the spectra for that gain
    '''    
    return DiffSpectra(sp1,gain*sp2,basebg)
#---------------------------------------------------------------------------
    
    
   
#----------------------------------------------------------------------------- 
#----------------------------------------------------------------------------- 
    
if __name__ == '__main__':

#-------------------------------------------------------------
 
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    for idx,arg in enumerate(sys.argv):
        if idx==1:
            num_arg=arg
        else:
            num_arg=None
            
    print 'number required = ', num_arg 
    
    if not num_arg.isdigit():
        num_arg=None           

    now=datetime.utcnow()  # choose UTC time
    datestr=str(now)
    print 'standard date format for the analysis :',datestr
    #  want the following format '2016-05-10T11:55:27.267'
    date_of_analysis=now.strftime('%Y-%m-%dT%H:%M:%S')
    print 'fits date format for the analysis : ',date_of_analysis

#------------------------------------------------------------------------


#---------------------------------------------------------------------------

    # input file : reduced and assembled images
    #--------------
    if num_arg==None:
        fileindex=1;
    else:
        fileindex=int(num_arg)
    
    object_name='HD158485_grat_'+str(fileindex)
    path='./HD158485_grat'
    basefilename='AssScImHD158485_grat_'+str(fileindex)+'.fits' # check master bias
    filename=os.path.join(path,basefilename)
    
    # output files
    #-------------
    #path
    spectrum_path=path+'_spectra'    # directory for spectra
      
    # basefilename
    basespectrumfilename='dataspec_'+object_name+'.fits'
    basespectrumfigfile='spec_'+object_name+'.pdf'
    
    # full output filename
    fullspectrumfilename=os.path.join(spectrum_path,basespectrumfilename)
    fullspecfigfilename=os.path.join(spectrum_path,basespectrumfigfile)
    
#-----------------------------------------------------------------------------
#   Selection Flags
#----------------------------------------------------------------------------
    RotationAngleOptimisation=False
    BackgroundSubtractionFlag=False
    CropLittleStarFlag=True

#--------------------------------------------------------------------------------

    #Get header
    hdulist = fits.open(filename)
    prim_hdr = hdulist[0].header
    exposure = prim_hdr['EXPOSURE']
    date_obs = prim_hdr['DATE-OBS']

    print date_obs
    print 'exposure = ',exposure,'seconds'

    # get image
    ccd_chan = ccdproc.CCDData.read(filename, hdu=0,unit='adu') 



#----------------------------------------------------------------------------------
# Show the original image
    fig, ax = plt.subplots(figsize=(8, 8))
    img=ax.imshow(ccd_chan,vmin=0,vmax=50.)
    plt.title(object_name)
    #plt.tight_layout()
    plt.colorbar(img)
    plt.grid(True)
    plt.show()
#------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# ## Background subtraction

# The user should decide if the background must be subtracted or not

    if BackgroundSubtractionFlag:
        bkg = Background2D(ccd_chan, (100, 100), filter_size=(3, 3),method='median')
        data=ccd_chan - bkg.background
        print('median background = ',bkg.background_median)
        print('median rms = ',bkg.background_rms_median)
    else:
        data=ccd_chan


    # normalisation to time
    data=data.divide(exposure)   # normalisation to one second exposure


    # This calculation of sky background has to be used to erase the central star

    skybackground=np.median(data)
    sigma_skybackground=np.std(data)


    print 'remaing sky background = {:2.3f} +/- {:2.3f}'.format(skybackground,sigma_skybackground)

#--------------------------------------------------------------------------------------
# ## 8.) Image rotation and spectrum region selection


    rotation_angle_test=np.linspace(67.0,67.5,100)
    NBTESTS=rotation_angle_test.shape[0]


# ### Optimisation to find the best rotation angle

    flux=np.zeros(NBTESTS)


    w=10   # spectrum width selection
           # ------------------------
    ws=100 # Central star width : must be that large especially for 5 seconds exposures
           #-------------------


    if RotationAngleOptimisation:
        for index,angle in np.ndenumerate(rotation_angle_test):
            rotated_image=ndimage.interpolation.rotate(data,angle)
            imax,jmax = np.unravel_index(rotated_image.argmax(), rotated_image.shape)
            region=rotated_image[imax-w: imax+w,480:4830]  # extract the region
            imax2,jmax2 = np.unravel_index(region.argmax(), region.shape)
            region[:,jmax2-ws:jmax2+ws]=skybackground  # suppress central star
            flux[index]=region.sum()
            print 'index = {} angle={:2.4f} flux = {:2.0f}'.format(index,angle,flux[index])


# In[882]:

    if RotationAngleOptimisation:
        dflux=flux-flux.max()
        fig, ax = plt.subplots(figsize=(8, 8))
        plt.plot(rotation_angle_test,dflux)
        plt.title(object_name)
        #plt.tight_layout()
        plt.colorbar(img)
        plt.grid(True)
        plt.show()
        


# In[883]:

    if not RotationAngleOptimisation:
        selected_angle=67.28787879


# In[884]:

    if RotationAngleOptimisation:
        selected_angle= rotation_angle_test[np.where(flux==flux.max())]


# In[885]:

    print 'selected angle = {} degrees'.format(selected_angle)


#---------------------------------------------------------------------------------
#  Rotation of the image


    rotated_image=ndimage.interpolation.rotate(data,selected_angle)



    #check the rotated image
    plt.figure(figsize=(8.,8.))
    plt.imshow(rotated_image,vmin=0,vmax=50)
    plt.show()


    wcheck=100  # defines the width
    imax,jmax = np.unravel_index(rotated_image.argmax(), rotated_image.shape)
    check_region=np.copy(rotated_image[imax-wcheck: imax+wcheck,510:4400])  # extract the region
 
    fig, ax = plt.subplots(figsize=(20, 3))
    ax.imshow(check_region,vmin=0,vmax=50)
    plt.show()    



    themaximum=rotated_image.max()



    imax,jmax = np.unravel_index(rotated_image.argmax(), rotated_image.shape)
    print imax,' ',jmax



    np.where(rotated_image==rotated_image.max())


# check the profiles

    profile1=np.sum(rotated_image[:,800:1000],axis=1)
    profile2=np.sum(rotated_image[:,1000:2000],axis=1)
    profile3=np.sum(rotated_image[:,3000:3800],axis=1)
    profile4=np.sum(rotated_image[:,4000:4600],axis=1)



    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col',sharey='row',figsize=(10.,8))
    ax1.semilogy(profile1)
    ax2.semilogy(profile2)
    ax3.semilogy(profile3)
    ax4.semilogy(profile4)
    #ax1.plot(profile1)
    #ax2.plot(profile2)
    #ax3.plot(profile3)
    #ax4.plot(profile4)
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax4.grid(True)


    the_imax1=np.where(profile1==profile1.max())
    the_imax2=np.where(profile2==profile2.max())
    the_imax3=np.where(profile3==profile3.max())
    the_imax4=np.where(profile4==profile4.max())

    print the_imax1,the_imax2,the_imax3,the_imax4,


    imax=int(np.median([the_imax1,the_imax2,the_imax3,the_imax4]))


    print imax # correct for the vertical position of the center
     #------------------------------------------------ 


    w=10  # This is the width we will use to get the spectrum
          #---------------------------------------------------



    # check the central region is OK
    # -------------------------------
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(rotated_image,vmin=0,vmax=50.)
    ax.plot([500, 4830], [imax-w, imax-w], color='y', linestyle='-', linewidth=2)
    ax.plot([500, 4830], [imax+w, imax+w], color='y', linestyle='-', linewidth=2)
    plt.title(object_name)
    #plt.tight_layout()
    plt.show()

#----------------------------------------------------------------------------------
# Big view of the spectrum region of reference
#----------------------------------------------

    # it also contains the region used to compute background

    SpectrumRegion=rotated_image[imax-100: imax+100,:]

    # Copy of the working region
    #-------------------------
    SpectrumRegionNew=np.copy(SpectrumRegion)   # must coy the original spectrum not to overwrite the original

    imax2,jmax2 = np.unravel_index(SpectrumRegionNew.argmax(), SpectrumRegionNew.shape)
    print imax2,' ',jmax2



    # put background where the central star is
    #-----------------------------------------
    SpectrumRegionNew[:,jmax2-ws: jmax2+ws]=skybackground # remove the big central star and replace by background


    # ### Extract the region of the spectrum  from the rotated image


    # show the original spectrum region
    #---------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SpectrumRegion,vmin=0,vmax=100.)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(True)
    plt.title(object_name,fontsize=30)
    plt.show()


    # check the vertical position of the center along the dispersion axis
    #---------------------------------------------------------------------
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col',figsize=(10.,10))
    ax1.plot(np.sum(SpectrumRegion[:,800:1000],axis=1))
    ax2.plot(np.sum(SpectrumRegion[:,1000:2000],axis=1))
    ax3.plot(np.sum(SpectrumRegion[:,3000:3800],axis=1))
    ax4.plot(np.sum(SpectrumRegion[:,4000:4600],axis=1))
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax4.grid(True)
    plt.show()


    # ### Show the spectrum region with the central star removed


    # show the spectrum cropped from central star
    #--------------------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SpectrumRegionNew,vmin=0,vmax=100.)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(True)
    plt.title(object_name,fontsize=30)
    plt.show()


    # idem but with a zoom
    #----------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SpectrumRegionNew[:,300:4500],vmin=0,vmax=50.)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(True)
    plt.title(object_name,fontsize=30)
    plt.show()

#-------------------------------------------------------------------------------
    # ### Crop little star on the right wing


    ymax3,xmax3=np.where(SpectrumRegionNew==SpectrumRegionNew.max())
    print ymax3,xmax3

    littlestarsize=10    # be carefull the extinction of the tiny star should not be that big


    # total size the of Spectrum region in X and Y
    #----------------------------------------------
    X_SRN=np.arange(SpectrumRegionNew.shape[1])
    Y_SRN=np.arange(SpectrumRegionNew.shape[0])


    # generate a grid (XV_SRN, YV_SRN are 2D)
    #----------------------------------------
    XV_SRN, YV_SRN = np.meshgrid(X_SRN, Y_SRN, sparse=False, indexing='ij')



    # Get the coordinates of the points inside the star disk
    #--------------------------------------------------------
    Xlsr,Ylsr=np.where( (XV_SRN-xmax3)**2+(YV_SRN-ymax3)**2<littlestarsize**2)
    selected_indexes=zip(Ylsr,Xlsr)



    # here the little star is erased
    #---------------------------------
    for ii,jj in selected_indexes:
        SpectrumRegionNew[ii,jj]=skybackground



    # check the liitle star is erased
    #---------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SpectrumRegionNew[:,3000:4500],vmin=0,vmax=100.)
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=10)
    plt.title(object_name,fontsize=30)
    plt.show()


    #### check the transverse profile
    part=SpectrumRegion[:,4000:4500]
    pro=np.sum(part,axis=1)


    # see la patte d'oiseau à la fin 
    plt.figure(figsize=(8.,8.))
    plt.plot(pro)
    plt.show()

#--------------------------------------------------------------------------------------------
# ## Optimize the width
# 
# 
# The width is varied from zero to the Spectrum Region width to study the "encircled energy" versus width


    widths=np.arange(1,100,1)
    NBWIDTHS=widths.shape[0]

    profile_list=[]
    integratedflux=np.zeros(NBWIDTHS)

    plt.figure(figsize=(8.,8.))
    plt.imshow(rotated_image,vmin=0,vmax=50) 
    plt.show()


    for index,thew in np.ndenumerate(widths):  
        TestRegion=np.copy(rotated_image[imax-thew: imax+thew,480:4830])
        imax2,jmax2 = np.unravel_index(TestRegion.argmax(), TestRegion.shape)
        TestRegion[:,jmax2-ws: jmax2+ws]=skybackground  # remove the star
        TestRegion=TestRegion-skybackground  # now must remove the sky background to all the pixels
        profile=np.sum(TestRegion,axis=1)
        integratedflux[index]=profile.sum()
        profile_list.append(profile)
    

    imax2,jmax2 = np.unravel_index(SpectrumRegion.argmax(), SpectrumRegion.shape)
    print imax2,' ',jmax2

    for index,thew in np.ndenumerate(widths):
    
        TestRegion=np.copy(SpectrumRegionNew[imax2-thew: imax2+thew,:])   # this spectrum has already the little star removed
        #imax2,jmax2 = np.unravel_index(TestRegion.argmax(), TestRegion.shape)
        #TestRegion[:,jmax2-ws: jmax2+ws]=skybackground  # remove the star
        TestRegion=TestRegion-skybackground  # now must remove the sky background to all the pixels
        profile=np.sum(TestRegion,axis=1)
        integratedflux[index]=profile.sum()
        profile_list.append(profile)
    
    integratedflux=integratedflux/integratedflux.max()

    plt.figure(figsize=(8.,8.))
    plt.plot(widths,integratedflux)
    plt.title('Spectrum transverse profile')
    plt.xlabel('number of pixels')
    plt.ylabel('fraction of signal')
    plt.show()
    

#-----------------------------------------------------------------------------------
# ### Select the final window width for the spectrum
    wsel=20  # should be 20, but set 10 due to the star


    # ## 9) Final Extraction of the spectrum from the image


# take again the Selected spectrum region but this time with the selected width wsel
#---------------------------------------------------------------------------------
    imax,jmax = np.unravel_index(rotated_image.argmax(), rotated_image.shape)
    print imax,' ',jmax
    SelectedSpectrumRegion=np.copy(rotated_image[imax-wsel: imax+wsel,:])



    # plot the final selected region for the spectrum
    #---------------------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SelectedSpectrumRegion,vmin=0,vmax=50.)
    plt.title(object_name)
    plt.show()


#--------------------------------------------------------------------------------
# Remove the central star in the new selected region
#-------------------------------------------------------
    imax2,jmax2 = np.unravel_index(SelectedSpectrumRegion.argmax(), SelectedSpectrumRegion.shape)
    print imax2,' ',jmax2

    SelectedSpectrumRegion[:,jmax2-ws: jmax2+ws]=skybackground


    # Check the central star has been erased
    # ---------------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SelectedSpectrumRegion,vmin=0,vmax=50.)
    plt.title(object_name)
    plt.tight_layout()
    plt.show()



#------------------------------------------------------------------------------
# Crop the littlestar of selected
#----------------------------------
    if CropLittleStarFlag:
        ymax4,xmax4=np.where(SelectedSpectrumRegion==SelectedSpectrumRegion.max())
        X_SRN=np.arange(SelectedSpectrumRegion.shape[1])
        Y_SRN=np.arange(SelectedSpectrumRegion.shape[0])
        XV_SRN, YV_SRN = np.meshgrid(X_SRN, Y_SRN, sparse=False, indexing='ij')
        Xlsr,Ylsr=np.where( (XV_SRN-xmax4)**2+(YV_SRN-ymax4)**2<littlestarsize**2)
        selected_indexes=zip(Ylsr,Xlsr)
        for ii,jj in selected_indexes:
            SelectedSpectrumRegion[ii,jj]=skybackground



    # Check the little star has been removed
    #----------------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(SelectedSpectrumRegion,vmin=0,vmax=50.)
    plt.title(object_name)
    plt.tight_layout()
    plt.show()


#---------------------------------------------------------------------------------
# ## Extract up and down bands to compute the background



    SelectedBackgroundUpRegion=np.copy(rotated_image[imax-wsel-100: imax+wsel-100,:])
    SelectedBackgroundDownRegion=np.copy(rotated_image[imax-wsel+100: imax+wsel+100,:])
    SelectedBackgroundMinRegion=np.where(np.less_equal(SelectedBackgroundUpRegion,SelectedBackgroundDownRegion),SelectedBackgroundUpRegion , SelectedBackgroundDownRegion)


    # plot the two background region
    # ---------------------------------------
    fig, (ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(20, 2))
    ax1.imshow(SelectedBackgroundUpRegion,vmin=0,vmax=10.)
    ax2.imshow(SelectedBackgroundDownRegion,vmin=0,vmax=10.)
    ax3.imshow(SelectedBackgroundMinRegion,vmin=0,vmax=10.)
    plt.tight_layout()
    plt.show()


    BackgUp=np.median(SelectedBackgroundUpRegion,axis=0)  # median along the vertical
    BackgDo=np.median(SelectedBackgroundDownRegion,axis=0) # median along the vertical
    BackgMi=np.median(SelectedBackgroundMinRegion,axis=0) # median along the vertical

    m1=np.median(BackgUp)
    m2=np.median(BackgDo)
    m0=np.median(BackgMi)

    BackgMi=BackgMi+ ((m1+m2)/2-m0)  # remove the bias induced by the min

    fig, ax = plt.subplots(figsize=(20, 5))
    ax.plot(BackgUp,color='m',label='background from up band')
    ax.plot(BackgDo,color='b',label='background of down band')
    ax.plot(BackgMi,color='r',label='background of min band')
    #plt.ylim(0,4.5)
    plt.legend()
    plt.show()

    FinalBackground=SelectedBackgroundMinRegion+((m1+m2)/2-m0)



    # plot the final background region
    # ---------------------------------------
    fig, ax1 = plt.subplots(1,1,figsize=(20, 2))
    ax1.imshow(FinalBackground,vmin=0,vmax=10.)
    plt.tight_layout()
    plt.show()



    BackgSpecFinal=np.median(FinalBackground,axis=0) # median along the vertical


    fig, ax = plt.subplots(figsize=(20, 5))
    ax.plot(BackgSpecFinal,color='r',label='Final Background')
    plt.legend()
    #plt.ylim(0,4.5)
    plt.show()

#--------------------------------------------------------------------------------
# Final spectrum

    FinalSpectrumRegion=SelectedSpectrumRegion-FinalBackground



    # Check the FinalSpectrumRegion
    #----------------------------------------
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.imshow(FinalSpectrumRegion,vmin=0,vmax=50.)
    plt.title(object_name)
    plt.tight_layout()
    plt.show()


#-------------------------------------------------------------------------------
# # Extract the 1D spectrum

    #spectrum=np.sum(SelectedSpectrumRegion,axis=0)
    spectrum=np.sum(FinalSpectrumRegion,axis=0)



    # ## 10.)  Plot the spectrum
    fig, ax = plt.subplots(figsize=(20, 8))
    #plt.semilogy(spectrum)
    plt.plot(spectrum)
    plt.ylim(-50.,2000.)
    plt.title(object_name)
    plt.tight_layout()
    #plt.xlim(1500.,2500.)
    plt.grid(True)


    # get the background of the 1D spectrum for chi2
    specbg=GetSpectrumBackground(spectrum,2000,3000,skybackground*wsel)
    
    plt.figure(figsize=(8.,8.))
    plt.plot(specbg)
    plt.show()
    

    bbg1=specbg[0:300]
    bbg2=specbg[700:1000]
    bbg3=[bbg1,bbg2]

    avspecbg=np.median(bbg3)
    rmsspecbg=specbg.std()
    print 'spectra noise = {} +/- {} '.format(avspecbg,rmsspecbg)


    # test the separation of the two spectra
    spec1,spec2=SeparateSpectra(spectrum,2500)


    fig, ax = plt.subplots(figsize=(20, 8))
    ax.plot(spec1)
    ax.plot(spec2)
    plt.title(object_name)
    plt.tight_layout()
    plt.ylim(-50.,2000.)
    plt.show()



    # Find the center between the two spectrum
    origins,thechi2=FindCenter(spectrum,2000,3000,rmsspecbg)

    plt.figure(figsize=(8.,8.))
    plt.plot(origins,thechi2)
    plt.ylim(0,200.)
    plt.show()



    indexmin=np.where(thechi2==thechi2.min())[0]
    theorigin=origins[indexmin]
    print indexmin[0],theorigin[0],thechi2.min()



    print 'the chi2 min',thechi2[515]


    # the two spectra are separated
    spec1,spec2=SeparateSpectra(spectrum,theorigin[0])



    fig, ax = plt.subplots(figsize=(20, 8))
    ax.plot(spec1,color='r',label='left spectrum')
    ax.plot(spec2,color='b',label='right spectrum')
    plt.title(object_name)
    plt.legend(loc='best')
    plt.ylim(-50.,1500.)
    plt.tight_layout()


    # ### Tiny amplitude correction
    # 
    # This correction must be applied, because I had not corrected for the differences in amplifier gain


    bbbg1=spec1[0:300]
    bbbg2=spec2[0:300]
    bbbg=0.5*(bbbg1.std()+bbbg2.std())
    print bbbg


    gains=np.linspace(0.5,1.5,1000)
    NBGAINS=gains.shape[0]
    chi2New=np.zeros(NBGAINS)

    # compute the chi2
    for idx,gain in np.ndenumerate(gains):
        chi2New[idx]=DiffAmplitudes(spec1,spec2,gain,bbbg)
    
    plt.figure(figsize=(8.,8.))
    plt.plot(gains,chi2New)
    plt.show()

    thegain=gains[np.where(chi2New==chi2New.min())]

    print 'thegain=',thegain


    # correct one of the spectrum
    spec2=spec2*thegain

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.plot(spec1,color='r',label='left spectrum')
    ax.plot(spec2,color='b',label='right spectrum')
    plt.title(object_name)
    plt.legend(loc='best')
    plt.ylim(-50.,1500.)
    plt.tight_layout()


    # Extract the final spectrum 
    N1=spec1.shape[0]
    N2=spec2.shape[0]
    N=np.min([N1,N2])
    spec1_croped=spec1[0:N]
    spec2_croped=spec2[0:N]
    specsum=0.5*(spec1_croped+spec2_croped)-avspecbg

    fig, ax = plt.subplots(figsize=(15, 8))
    ax.plot(specsum)


    # O2 Frauwnoffer A 759.370
    ax.plot([1315, 1315], [-200,2000], color='r', linestyle='-', linewidth=2)
    ax.plot([1400, 1400], [-200,2000], color='r', linestyle='-', linewidth=2)

    # H-alpha : Hα     656.281 nm
    ax.plot([1120, 1120], [-200,2000], color='m', linestyle='-', linewidth=2)
    ax.plot([1200, 1200], [-200,2000], color='m', linestyle='-', linewidth=2)

    # H-beta : Hβ     486.134 nm
    ax.plot([830, 830], [-200,2000], color='g', linestyle='-', linewidth=2)
    ax.plot([900, 900], [-200,2000], color='g', linestyle='-', linewidth=2)

    # H-gamma : Hγ     434.047 nm
    ax.plot([750, 750], [-200,2000], color='c', linestyle='-', linewidth=2)
    ax.plot([790, 790], [-200,2000], color='c', linestyle='-', linewidth=2)


    # H-delta : Hδ     410.175 nm
    ax.plot([715, 715], [-200,2000], color='k', linestyle='-', linewidth=2)
    ax.plot([740, 740], [-200,2000], color='k', linestyle='-', linewidth=2)

    # H-epsilon : H epsilon 397,0 nm
    ax.plot([695, 695], [-200,2000], color='y', linestyle='-', linewidth=2)
    ax.plot([710, 710], [-200,2000], color='y', linestyle='-', linewidth=2)

    plt.ylim(0.,2000.)
    plt.title(object_name)
    plt.tight_layout()


#-------------------------------------------------------------------------------
# Wavelength calibration
#--------------------------------------------------------------------------------

    peak_O2=np.min(specsum[1315:1400])
    indexes_O2=np.where(specsum==peak_O2)
    print indexes_O2, peak_O2


    peak_Halpha=np.min(specsum[1120:1200])
    indexes_Halpha=np.where(specsum==peak_Halpha)
    print indexes_Halpha, peak_Halpha


    peak_Hbeta=np.min(specsum[830:900])
    indexes_Hbeta=np.where(specsum==peak_Hbeta)
    print indexes_Hbeta, peak_Hbeta



    peak_Hgamma=np.min(specsum[750:790])
    indexes_Hgamma=np.where(specsum==peak_Hgamma)
    print indexes_Hgamma, peak_Hgamma


    peak_Hdelta=np.min(specsum[715:740])
    indexes_Hdelta=np.where(specsum==peak_Hdelta)
    print indexes_Hdelta, peak_Hdelta


    peak_Hepsilon=np.min(specsum[695:710])
    indexes_Hepsilon=np.where(specsum==peak_Hepsilon)
    print indexes_Hepsilon, peak_Hepsilon


    # X-axis of wavelength calibration line
    pixel_axis=np.array([indexes_Hepsilon[0] ,indexes_Hdelta[0], indexes_Hgamma[0], indexes_Hbeta[0],indexes_Halpha[0], indexes_O2[0]])

    # Y-axis of wavelength calibration line
    wavelength_axis=np.array([397.0, 410.175, 434.047,486.134, 656.281, 759.370])



    fig, ax = plt.subplots(figsize=(10, 10))
    plt.plot(pixel_axis,wavelength_axis,'o-')
    plt.plot([0,1351],[0,759.370],'k-')
    plt.title('Ronchi Grating Calibration line')
    plt.xlabel('number of pixels')
    plt.ylabel('wavelength in nm')
    plt.ylim(0,1000.)
    plt.xlim(0,1500.)
#plt.annotate('$O_2$', xy=(1200,763.370), xytext=(1000,763.37),fontsize=30,
#            arrowprops=dict(facecolor='black', shrink=0.02))
#plt.annotate('$H_{alpha}$', xy=(1200,658.281), xytext=(1000,658.281),fontsize=30,
#            arrowprops=dict(facecolor='black', shrink=0.02))
#plt.annotate('$H_{beta}$', xy=(1100,486.281), xytext=(900,486.281),fontsize=30,
#           arrowprops=dict(facecolor='black', shrink=0.02))
#plt.annotate('$H_{gamma}$', xy=(800,434.281), xytext=(600,385.281),fontsize=30,
#           arrowprops=dict(facecolor='black', shrink=0.02))
#plt.annotate('$H_{delta}$', xy=(600,410.500), xytext=(400,410.281),fontsize=30,
#            arrowprops=dict(facecolor='black', shrink=0.02))







    pixel_to_wavelength_spl = UnivariateSpline(pixel_axis,wavelength_axis)
    wavelength_to_pixel_spl = UnivariateSpline(wavelength_axis,pixel_axis)



    pixel_to_wavelength_spl(indexes_O2[0])
    index_O2=np.where(specsum==peak_O2)


    # fit for calibration line

    #np.polyfit(pixel_axis[:,0],wavelength_axis,1,full=True)
    np.polyfit(pixel_axis[:,0],wavelength_axis,1)


    # get the wavelength array of same size as the spectrum
    specsum_indexes=np.arange(specsum.shape[0])   
    specsum_wavelength=pixel_to_wavelength_spl(specsum_indexes)


    # Plot the calibrated experimental spectrum
    plt.figure(figsize=(15.,8.))
    plt.plot(specsum_wavelength,specsum)
    plt.title('Spectrum reconstructed with Monocam and Ronchi Grating',fontsize=30)
    plt.xlabel('$\lambda$ (nm)',fontsize=20)
    plt.ylabel(' flux in arbitrary units ',fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.ylim(-100.,1500.)
    plt.annotate('$O_2$', xy=(763.370,900), xytext=(763.37, 1500),fontsize=30,
            arrowprops=dict(facecolor='black', shrink=0.02))
    plt.annotate('$H_{alpha}$', xy=(658.281,1200), xytext=(658.281, 1500),fontsize=30,
            arrowprops=dict(facecolor='black', shrink=0.02))
    plt.annotate('$H_{beta}$', xy=(486.281,1300), xytext=(486.281,1500),fontsize=30,
           arrowprops=dict(facecolor='black', shrink=0.02))
    plt.annotate('$H_{gamma}$', xy=(434.281,1000), xytext=(385.281, 1500),fontsize=30,
           arrowprops=dict(facecolor='black', shrink=0.02))
    plt.annotate('$H_{delta}$', xy=(410.500,600), xytext=(410.281,200),fontsize=30,
            arrowprops=dict(facecolor='black', shrink=0.02))
#wavelength_axis=np.array([397.0, 410.175, 434.047,486.134, 656.281, 759.370])
    plt.ylim(0.,2000.)
    plt.xlim(0.,1000.)

#-----------------------------------------------------------------------------------
#SED
#----------------------------------------------------------------------------------

    path_sed='../SED'
    obj_name='hd158485'
    airmass='1.1'
    night_name='20160509'
    tablefitsfile_sed='SEDPred_'+obj_name+'_'+night_name+'_'+str(fileindex)+'.fits'


    # file of SED
    fullfilename_sed=os.path.join(path_sed,tablefitsfile_sed)

    # get the content of the file
    hdulist=fits.open(fullfilename_sed)
    hdulist.info()

    table_data=hdulist[1].data
    table_data.columns  # shows the columns names of the table

    # retrieve the table inside SED file

    wavelength_sed=table_data.field('WAVELENGTH')
    flux_sed=table_data.field('SEDcalspec')
    flux_sedccd=table_data.field('SEDxQE')
    flux_sedccdatm=table_data.field('SEDxQExATM')
    flux_sedccdatmopt=table_data.field('SEDxQExATMxTopt')


    # rename
    SED1=flux_sed             # SED of the star
    SED2=flux_sedccd          # SED of star multipled by CCD efficiency 
    SED3=flux_sedccdatm       # SED of star multipled by CCD efficiency and atmosphere transparency
    SED4=flux_sedccdatmopt    # SED of star multipled by CCD efficiency,atmosphere transparency,
                              # optic throuput



    SED5=SED4

    basefigfile=""

    plt.figure(figsize=(15,8))
    plt.plot(wavelength_sed,SED1,label='calspec SED',color='k')
    plt.plot(wavelength_sed,SED2,label='with CCD-QE',color='b')
    plt.plot(wavelength_sed,SED3,label='with CCD-QE and atmosphere',color='r')
    plt.plot(wavelength_sed,SED4,label='with CCD-QE,Optics and atmosphere',color='g')
    plt.xlim(0,1200.)
    plt.ylim(0,SED1.max()*1.2)
    plt.title('prediction for SED observation for hd158485 (calspec)',fontsize=40)
    plt.xlabel('wavelength $\lambda$ (nm)',fontsize=30)
    plt.ylabel('Flux x Transmission ',fontsize=30)
    plt.tick_params(axis='x', labelsize=30)
    plt.tick_params(axis='y', labelsize=30)
    plt.legend(fontsize=30)
    plt.show()



    plt.figure(figsize=(15.,8.))
    plt.plot(specsum_wavelength,specsum,label='data',color='b')
    plt.plot(wavelength_sed,SED5*3e14,label='predicted SED(ccd,atm,opt)',color='r')
    thetitle='Spectrum reconstructed with Monocam and Ronchi Grating for '+object_name
    plt.title(thetitle,fontsize=30)
    plt.xlabel('$\lambda$ (nm)',fontsize=20)
    plt.ylabel(' flux in arbitrary units ',fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    plt.legend(loc='best')
    plt.xlim(0.,1000.)
    
    plt.savefig(fullspecfigfilename)
    plt.show()
    
    

##-------------------------------------------------------------------------------
# ## save spectrum in a file 
##-----------------------------------------------------------------------------

    #path
    #spectrum_path='./HD158485_grat'+'_spectra'

    # basefilename
    #basespectrumfilename='dataspec_'+object_name+'.fits'

    # full output filename
    #fullspectrumfilename=os.path.join(spectrum_path,basespectrumfilename)


    # keep the original header
    primheader=prim_hdr


    # add the object name
    primheader['OBJNAME']=object_name


    # define the primary fits unit
    primhdu=fits.PrimaryHDU(header=primheader)


    # define the columns of the fits table
    col1 = fits.Column(name='WAVELENGTH', format='E', array=specsum_wavelength)
    col2 = fits.Column(name='SPECTRUMDATA', format='E', array=specsum)


    # define the fits Extension unit for the spectrum data
    cols = fits.ColDefs([col1, col2])     # definition of the column
    tbhdu = fits.BinTableHDU.from_columns(cols)     # new binary table HDU



    # save data and SED prediction
    thdulist = fits.HDUList([primhdu, tbhdu,hdulist[1]])
    thdulist.writeto(fullspectrumfilename,clobber=True)

#--------------------------------------------------------------------------------
    print '=================================================================='
    print 'END'
    print '=================================================================='
 #-----------------------------------------------------------------------------   