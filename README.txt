Arthur Argall - E6891 Project

The script proj_6891_2 conducts a variety of data import and preprocessing methods. It uses the custom functions
tcga_filerenamer, tcga_arraymaker_2, and cel_extractor. The function FileRename from Jan Simon is also necessary,
and can be downloaded here: http://www.mathworks.com/matlabcentral/fileexchange/29569-filerename 

FileRename needs to be compiled, or a compiled mex file needs to be added to the MATLAB directory. The mex 
file can be downloaded here: http://www.n-simon.de/mex/

Datasets come from TCGA and GEO. I have uploaded the TCGA datasets as they are difficult to obtain. The GEO
datasets can be downloaded from the following links:

Breast: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034
Colon: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14333
Ovarian: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9891

At these links download the raw data (bottom of the page). These files are very large, and when unpacked by the
script proj_6891_2 the unzipped files take up ~16 GB of space.

Both tcga_filerenamer and cel_extractor require directory addresses as inputs. They are currently set to my own
computer's locations, so all dr variables in the script will need to be changed to where you download the data.

Sections 4 and 5 of the script require the MATLAB bioinformatics toolbox and an internet connection to download
the data from GEO. Section 4 will probably be phased out and section 5 may be as well, so if you do not have
the bioinformatics toolbox you will hopefully still be able to run the code.

As a general matter the script is currently quite slow and memory intensive. If your computer is crashing try 
commenting out section 4 and working with only one dataset at a time (ie only run tcga_arraymaker_2 once in section
3).

Finally, section 6/cel_extractor will only unzip files to a given directory a single time. If you need to unzip
the files again, choose a new output directory or delete all the files within the previous output directory.