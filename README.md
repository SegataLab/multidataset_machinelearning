# Scope and Requirements #

## This repo provides a series of script for Data-Driven Whole-Genome-Shotgun (WGS) Metagenomics, in particular:

* 1 program allowing the download of a dataset from NCBI.
* 2 programs for speed up the boring procedure of writing in Bash a whole Machine-Learning (ML) meta-analyses based on several WGS datasets (all or part of them, maybe, just downloaded with the above named program).  
* 3 programs producing figures relative to the different aspect of the ML analyses conducted with the above program (the indeed become quite detailed).
* 4 programs devoted to the basic part of a standard WGS analyses, including the generation of figures such as: any type of ordination plot (literally), beta-diversity box-plots, alpha diversity boxplots.
* 5 program plotting the correlation between an index and the distance from the centroid of a class in the data.

### ALL the listed programs can be quite tricky to run if never done before, 'cos a detailed documentaion is still lacking: for any necessity fell free to write at paolomanghi1974@gmail.com.
## generation of figures relative to the Machine-Learning analysis in the paper: "Metagenomic analysis of colorectal cancer datasets identifies cross-cohort microbial diagnostic signatures and a link with choline degradation." 

## Getting Started

To clone the present repo:

```
hg clone https://<USER>@bitbucket.org/CibioCM/multidataset_machinelearning
```

The analysis listed here are mainly dovoted to generate the figures of the paper. Can be applied in a similar fashions to other series of
datasets/other classification problems.

###### NOTE: they are THEY ARE BASED ON TWO OTHER PUBLIC REPOS: cmdpy and metaml.

To clone cmdpy and metaml repos:

```
hg clone https://<USER>@bitbucket.org/CibioCM/metaml
hg clone https://<USER>@bitbucket.org/CibioCM/cmdpy

```

cmdpy repo contains the tool cmdpy_dataset.py which is used to merge quantitative profiles with metadata
metaml repo contains the tool classification_thomas-manghi.py which is needed to run the main analysis

###### NOTE: these figures have been create from a server in which I'm also maintaining the curatedMetagenomicDataset package
###### that means that is conceived as an architecture with the following tree-shape:

- /base_folder
- /base_folder/dataset_name(e.g. FengQ_2015)
- /base_folder/dataset_name/dataset_name_metadata.tsv
- /base_folder/dataset_name/metaphlan2
- /base_folder/dataset_name/metaphlan2/sample_#1_name
- /base_folder/dataset_name/metaphlan2/sample_#2_name
- /base_folder/dataset_name/metaphlan2/sample_#1_name/sample_#1_name_profile.tsv
- /base_folder/dataset_name/metaphlan2/sample_#1_name/sample_#2_name_profile.tsv

##### once the architecture is set (that means also
##### having dataset-based folders, metadata table and metaphlan2 profiles)
##### the first step is to proceed to the dataset generation (ML == starting datasets)
##### THEREFORE, in order to generate the datasets: do

    ```
	python run.py crc --define study_condition:CRC:control \
		-ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 \
		-db metaphlan \
		-do cross_study \
			-g0 nt:500 -g1 nsl:5 -g2 c:entropy
	```

###### if run without errors, you should be able to see in the folder you run it to files.sh:

    ```
	ls
	
	datasetcaller_metaphlan.sh
	transfer_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
	lodo_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
	```

###### these .sh are wrappers for all the analysis needed
###### IF YOU INSTALLED THE CMDPY AND METAML REPOS YOU SHOULD BE ABLE TO RUN THESE COMMANDS IN THE ORDER LISTED

### Once the analysis have run, in order to generate the figures you might run:

	* figure 2 panel a:
	
	** This analysis consist of a cross prediction matrix consisting of 
	** n. datasets * n. datasets tests + n. datasets Leave-One-Dataset-Out predictions

	```
    	
	
	
	```
	* Figure 3 panel a:
	** This analysis consist in a feature-ranking of 7 cross-validation
	*** Therefore the step n. 1 is to perform these cross-validation
	*** In order to speed up the analysis, I built a code for naming files:
	
	** 
	
	
	```
	python ../../../multidataset_machinelearning/run.py crc \
			--define study_condition:CRC:control \
			--datasets \
					FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 \
			-al rf \
			-do heat_map 
			-db metaphlan 
			-g0 c:entropy -g1 nt:1000 -g2 nsl:5 -nif 5 
			--path Fig_Three/
	
	```



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