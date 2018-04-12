# Scope and Requirements #
 
* multidataset_machinlearning is a collection of Python programs aimed to speed up the most essential steps in a WGS metagenomic analyses.
* it relies upone several commonly used Python libraries: among these numpy, scipy, pandas, scikit-learn, matplotlib, seaborn & biopython are not part of the original native libraries.
* it is not meant to be used as a standalone (see paragraphs below).

## There are essentially 3 goals pursued by the script in the collection: ##
### Downloading and Curating public Whole-Genome-Shotgun cohorts ###

* ncbi_downloader.py & ncbi_downloader_dev.py are the programs for executing the physical download of public datasets from NCBI: these script can be applid to different data types, such as AMPLICON sequencing and Viromes.
* cure.py is a script containing many utility functionalities aimed to manually curate clinical metadata. Examples of the utility of this curation are the metadata table of the CuratedMetagenomicDataset R package.

### Performing explorative analyses on a meta-cohort or a meta-dataset ###

* cmdpy_dataset is a tool for collecting a theoretically unlimited number of datasets: it assumes that for any dataset involved in the meta-analyses, community profiles are available, so as corresponding metadata tables. The profiles for which the tool has been conceived were Metaphlan2 and Humann2 profiles, along the lines of the curatedMetagenomicDataset R package. It has been then applied to other types of profiling. 
* cmdpy_betadiversity.py: the data collected by cmdpy_dataset.py can be send in STDIN to cmdpy_betadiversity.py (see below).

#### Cmdpy_betadiversity.py ####

* this tool computes beta and alpha diversity, and provides graphical and statistical representation of the required result, including full-optional ordination plots (MDS, PCA, NMDS).

### A comprehensive machine learning meta-analyses ###

* 
* 
* 

* Downloading and curating public 
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact