# MonoCamMay2016
#----------------

- created May 26th 2016
- Updated June 8th 2016
- author Sylvie Dagoret-Campagne


## requirements
- numpy 
- matplotlib
- astropy
- ccdproc

(may be pandas)


And DS9
- pyds9


I use python for Anaconda, python version 2.7. For me:

- source activate pyastrophys



## Directories
------------

### Directory 1 : Exploration
---------------
Basic ipython notebook  to have a quick look to the data.
Not a detailed analysis

#### Look at single files

- ExploreBias.ipynb		
- ExploreDark.ipynb		
- ExploreFlats.ipynb		
- ExploreImages.ipynb		
- ExploreProjector.ipynb	
- ExploreSkyFlats.ipynb


#### scan over files of the same type

- ScanOverBiasFiles.ipynb
- ScanOverDarkFiles.ipynb
- ScanOverFlats.ipynb
- ScanOverSkyFlats.ipynb


#### tools

- MakeFileList.ipynb
- SaveFitFile.ipynb
- ShowFileInDS9.ipynb
- TestFitHisto_v0.ipynb


### Directory 2 : MyReduction_2016_05_09
 
#### very important tasks for standard CCD reduction:

- 2.a) MONOCAMCCDRed.ipynb --> Produce the Master bias and the Master Darks
- 2.b) MONOCAMCCDRedWithSkyFlats.ipynb --> Produce the Master Flat with Flat Sky

- 2.c) Work with the scientific images without assembling them

-- MONOCAMCCDRedImages_HD158485_grat.ipynb
-- MONOCAMCCDRedImages_HD159222_grat.ipynb
-- MONOCAMCCDRedImages_HD163466_grat.ipynb


- 2.d) Does the images assembly into a single 4K x 4K image: important to reconstruct the spectrum

-- ImageAssembler_HD158485_grat.ipynb
-- ImageAssembler_HD159222_grat.ipynb
-- ImageAssembler_HD163466_grat.ipynb


- 2.e) To study the sky background and check if we understand the amplifier gain

-- SkyBackgroundAndGain_HD158485_grat.ipynb
-- SkyBackgroundAndGain_HD159222_grat.ipynb
-- SkyBackgroundAndGain_HD163466_grat.ipynb



#### Python Library:

provide services and tools to ipython notebooks:

- libMonocamBaseImages.py

#### Example to produce a Master Bias and Master Dark


- Build_MasterBias.ipynb		
- Build_MasterDark.ipynb		

#### Utility tools

- ShowFileInDS9.ipynb
- ViewSingleFitsFile.ipynb
- WorkWithDS9.ipynb
- grab_files.ipynb


## To get the whole Package
---------------------------

https://github.com/sylvielsstfr/MonoCamMay2016.git

## To make a push do:
------------------------

git remote set-url origin git@github.com:sylvielsstfr/MonoCamMay2016.git
git push origin master
