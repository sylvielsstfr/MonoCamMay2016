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

### Exploration
---------------
Basic ipython notebook  to have a look to the data

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


### MyReduction_2016_05_09
 
#### very important tasks for standard CCD reduction:

- MONOCAMCCDRed.ipynb --> Produce the Master bias and the Master Darks
- MONOCAMCCDRedWithSkyFlats.ipynb __> Produce the Master Flat with Flat Sky


#### Python Library

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
