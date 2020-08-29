Project done by Svetlana Semenova.
Mentored by Dinara Usmanova from Columbia University, Department of Systems Biology.
Additional support provided by Angelique Bosse from Montgomery Blair High School.

analysis_raw_data.py asks for the location of a folder (from home directory) where raw expression data
and evolution rates reside. Then allows for saving and creating linear regression models, along with
other methods of statistical analysis.

Analysis has two folders: Partials and Images. Partials were only necessary for two datasets, so the 
folder is sparse. Images contains a folder for each dataset, where images were made for R^2 values, 
similarity between conditions values, and more. Original datasets were not used, only what 
analysis_raw_data.py outputed.

Data contains three folders: Escherichia_coli, Escherichia_coli_REL606, and Mycobacterium_tuberculosis.
The Mtb folder contains four distinct datasets. Within each dataset folder, the citations is in a text
file. All evolution rates were calculated using data from Ensembl Bacteria and using the Biopython
library.

Further documentation in progress.