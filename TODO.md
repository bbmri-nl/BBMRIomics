Here are a few things that were not high on the priority list but
could be beneficial for the *BBMRIomics*-package:

### Adding unit-test ###

R has several ways to incorporate unit-tests to packages which can
even be added to GitHub CI hook on pushing. 

We could think of checking existence of data files and the
accessibility of the database urls.

### Adding new updated databasets for Freeze2 ###

The metadata for the freeze2 samples has been updated since the original 
SummarizedExperiments were made for the RNAseq counts and methylation beta 
and M values. New SummarizedExperiments should be contructed with the updated
data.